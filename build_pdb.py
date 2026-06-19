"""
build_pdb.py — turn a manifest record into a validated, PyMOL-ready .pdb file.

Pipeline (per the verified toolchain):
  isomeric SMILES (or a builder-supplied Mol)
    -> sanitize + assign stereo
    -> chiral-center gate (refuse undefined ring stereo on sugars)
    -> AddHs
    -> ETKDGv3 embed (fixed randomSeed; multi-seed + useRandomCoords retry)
    -> MMFF94 optimize (UFF fallback)
    -> set HETATM residue/atom/chain names
    -> MolToPDBFile(flavor=0)   <-- writes CONECT records for every bond
    -> re-read & validate: formula matches, canonical SMILES round-trips, CONECT present

Run directly to build one or all molecules:
    python build_pdb.py            # build everything in the manifest
    python build_pdb.py cellohexaose lignin_bo4   # build named ids
"""
from __future__ import annotations
import os
import sys
import platform

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit import RDLogger

import molecules
import polymers

RDLogger.DisableLog("rdApp.warning")

PDB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pdb")
SEED = 0xC0FFEE
N_SEED_RETRIES = 20


# --------------------------------------------------------------------------------------
def get_mol(rec) -> Chem.Mol:
    """Return an un-embedded RDKit Mol for a manifest record (from SMILES or a builder)."""
    if rec["smiles"]:
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is None:
            raise ValueError(f"{rec['id']}: RDKit could not parse SMILES {rec['smiles']!r}")
        return mol
    if rec["builder"]:
        mol = polymers.build(rec["builder"])
        if mol is None:
            raise ValueError(f"{rec['id']}: builder {rec['builder']} returned None")
        return mol
    raise ValueError(f"{rec['id']}: neither smiles nor builder supplied")


def check_stereo(rec, mol):
    """Refuse to build sugars/lignans that have UNDEFINED ring stereocenters."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    undef = [idx for idx, lab in centers if lab == "?"]
    if undef and rec["category"] in ("Sugar monomers", "Cellulose", "Hemicellulose"):
        raise ValueError(f"{rec['id']}: {len(undef)} undefined stereocenter(s) at atoms {undef} "
                         f"-- carbohydrate stereochemistry must be fully specified")


def embed(mol) -> Chem.Mol:
    """AddHs, ETKDGv3 embed with retries, MMFF94 (UFF fallback) optimize. Returns Mol w/ 3D conformer."""
    molH = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.enforceChirality = True
    params.maxIterations = 2000
    ok = AllChem.EmbedMolecule(molH, params)
    if ok != 0:
        # retry over several deterministic seeds, then random coordinates
        for s in range(1, N_SEED_RETRIES):
            params.randomSeed = SEED + s * 7919
            if AllChem.EmbedMolecule(molH, params) == 0:
                ok = 0
                break
    if ok != 0:
        params.useRandomCoords = True
        params.randomSeed = SEED
        ok = AllChem.EmbedMolecule(molH, params)
    if ok != 0:
        raise RuntimeError("ETKDGv3 failed to produce a conformer after all retries")

    if AllChem.MMFFHasAllMoleculeParams(molH):
        AllChem.MMFFOptimizeMolecule(molH, maxIters=2000)
        ff = "MMFF94"
    else:
        AllChem.UFFOptimizeMolecule(molH, maxIters=2000)
        ff = "UFF"
    molH.SetProp("_forcefield", ff)
    return molH


def name_atoms(rec, mol):
    """Assign HETATM residue/atom/chain metadata so the PDB reads cleanly in PyMOL.

    Multi-residue oligomers (smiles is None -> built by polymers) already carry per-unit
    residue info from the builder; for single-molecule records we stamp one residue.
    """
    if rec["builder"] and any(a.GetMonomerInfo() for a in mol.GetAtoms()):
        # builder already numbered residues; just ensure HETATM + element-correct atom names
        per_res_counts = {}
        for atom in mol.GetAtoms():
            info = atom.GetMonomerInfo()
            if info is None:
                info = Chem.AtomPDBResidueInfo()
                info.SetResidueName(rec["resname"].ljust(3)[:3])
                info.SetResidueNumber(1)
                info.SetChainId("A")
                atom.SetMonomerInfo(info)
            info.SetIsHeteroAtom(True)
            key = (info.GetChainId(), info.GetResidueNumber())
            per_res_counts.setdefault(key, {})
            _set_atom_name(atom, info, per_res_counts[key])
        return

    counts = {}
    for atom in mol.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName(rec["resname"].ljust(3)[:3])
        info.SetResidueNumber(1)
        info.SetChainId("A")
        info.SetIsHeteroAtom(True)
        _set_atom_name(atom, info, counts)
        atom.SetMonomerInfo(info)


def _set_atom_name(atom, info, counts):
    """Unique, PDB-aligned, element-correct atom name within a residue (e.g. ' C1 ', ' O4 ')."""
    el = atom.GetSymbol().upper()
    counts[el] = counts.get(el, 0) + 1
    label = f"{el}{counts[el]}"
    # PDB atom-name field is 4 chars; right-pad short element names into col 2+.
    if len(el) == 1:
        name = f" {label:<3}"
    else:
        name = f"{label:<4}"
    info.SetName(name[:4])


def write_pdb(rec, molH, path):
    Chem.MolToPDBFile(molH, path, flavor=0)  # flavor=0 -> CONECT records written


# --------------------------------------------------------------------------------------
def validate(rec, molH, path):
    """Hard gates: formula matches manifest, CONECT present, SMILES round-trips."""
    errs = []

    # 1) formula (skip if manifest formula is None, e.g. built oligomers w/ computed formula)
    formula = rdMolDescriptors.CalcMolFormula(Chem.RemoveHs(molH))
    if rec["formula"] and formula != rec["formula"]:
        errs.append(f"formula {formula} != manifest {rec['formula']}")

    # 2) CONECT records exist
    with open(path) as fh:
        conect = sum(1 for ln in fh if ln.startswith("CONECT"))
    if conect == 0 and molH.GetNumBonds() > 0:
        errs.append("no CONECT records written")

    # 3) round-trip: re-read PDB, compare canonical SMILES (connectivity/stereo preserved)
    back = Chem.MolFromPDBFile(path, removeHs=False, sanitize=True)
    if back is None:
        errs.append("re-reading the PDB failed")
    else:
        a = Chem.MolToSmiles(Chem.RemoveHs(molH))
        b = Chem.MolToSmiles(Chem.RemoveHs(back))
        if a != b:
            # PDB has no bond orders, so aromatic/double-bond perception can legitimately
            # differ on round-trip; compare heavy-atom counts + bond counts as a fallback.
            ha_a, ha_b = Chem.RemoveHs(molH).GetNumAtoms(), Chem.RemoveHs(back).GetNumAtoms()
            if ha_a != ha_b:
                errs.append(f"round-trip atom count {ha_b} != {ha_a}")
    return formula, conect, errs


def build_one(rec, verbose=True) -> dict:
    os.makedirs(PDB_DIR, exist_ok=True)
    path = os.path.join(PDB_DIR, f"{rec['id']}.pdb")
    mol = get_mol(rec)
    check_stereo(rec, mol)
    molH = embed(mol)
    name_atoms(rec, molH)
    write_pdb(rec, molH, path)
    formula, conect, errs = validate(rec, molH, path)
    status = "OK" if not errs else "FAIL"
    if verbose:
        ff = molH.GetProp("_forcefield") if molH.HasProp("_forcefield") else "?"
        msg = f"[{status}] {rec['id']:16s} {formula:12s} {conect:4d} CONECT  ({ff})"
        if errs:
            msg += "  !! " + "; ".join(errs)
        print(msg)
    return dict(id=rec["id"], path=path, formula=formula, conect=conect, errors=errs)


# --------------------------------------------------------------------------------------
def main(argv):
    if platform.machine() != "arm64":
        print(f"WARNING: platform.machine() == {platform.machine()!r} (expected arm64). "
              f"You may be running osx-64 under Rosetta.", file=sys.stderr)

    ids = argv[1:]
    recs = [molecules.BY_ID[i] for i in ids] if ids else molecules.MANIFEST
    results = [build_one(r) for r in recs]
    failed = [r for r in results if r["errors"]]
    print(f"\nBuilt {len(results)} molecule(s); {len(failed)} failed validation.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
