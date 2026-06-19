"""
polymers.py — assemble cellulose / softwood-hemicellulose / lignin oligomers as RDKit Mols.

STRATEGY (stereo-safe SMILES ring-bond joining)
-----------------------------------------------
Each monomer is emitted as a SMILES *fragment* that is byte-for-byte the verified monomer
from molecules.MONOMER_SMILES, except that:
  * a glycosidic-donor anomeric carbon has its -OH replaced IN PLACE by a ring-closure label
    (e.g. `[C@@H](O)` -> `[C@@H]%10`), and
  * a glycosidic-acceptor hydroxyl gets the SAME ring-closure label appended (`O` -> `O%10`).
Fragments are concatenated with '.' ; the shared labels stitch them into ONE connected
molecule. Because we only swap a substituent token for a ring-bond *at the same position*,
no chiral atom's neighbor ORDER changes, so every stereocenter is preserved exactly.
A correctness gate (`_assert_templates`) confirms each all-free fragment canonicalizes to the
verified monomer, and `build("cellulose:2")` is checked against known beta-cellobiose elsewhere.

Real bonds => RDKit writes CONECT records for every linkage (what PyMOL needs).

dispatch tokens (the manifest's `builder` field):
    cellulose:N      N-mer beta-1,4 D-glucose ribbon (N in 2..8)
    ggm:compact      galactoglucomannan, 6-unit backbone + 1 Gal branch
    ggm:extended     10-unit backbone + 2 Gal branches + O-acetyls
    agx:compact      arabinoglucuronoxylan, 5-xylose backbone + MeGlcA + Araf
    agx:extended     8-xylose backbone + 2 MeGlcA + 2 Araf
    lignin:trimer    defined guaiacyl beta-O-4 / beta-5 trimer
"""
from __future__ import annotations
from rdkit import Chem

import molecules

# Glycosidic ring-bond labels start at 10 so they never collide with monomers' internal `1`.
_L0 = 10


def _lab(n):
    return f"%{n}"


# --------------------------------------------------------------------------------------
# Monomer fragment emitters. Each free slot defaults reproduce the verified monomer SMILES.
# A slot set to an int => that position carries a glycosidic ring-bond with that label.
# --------------------------------------------------------------------------------------
def _hexose(kind, c1=None, o4=None, o6=None, o2=None, o3=None, ac2=False, ac3=False):
    """beta-D-glucopyranose ('glc') / beta-D-mannopyranose ('man') fragment.
    Templates differ only at the C2 stereocenter."""
    o6s = f"O{_lab(o6)}" if o6 is not None else "O"
    c1s = _lab(c1) if c1 is not None else "(O)"
    o2s = "(OC(C)=O)" if ac2 else (f"(O{_lab(o2)})" if o2 is not None else "(O)")
    o3s = "(OC(C)=O)" if ac3 else (f"(O{_lab(o3)})" if o3 is not None else "(O)")
    o4s = f"O{_lab(o4)}" if o4 is not None else "O"
    # glucose: C2 is [C@H]; mannose (C2 epimer): C2 is [C@@H]
    c2 = "[C@H]" if kind == "glc" else "[C@@H]"
    return f"{o6s}C[C@H]1O[C@@H]{c1s}{c2}{o2s}[C@@H]{o3s}[C@@H]1{o4s}"


def _gal_donor(c1):
    """alpha-D-galactopyranosyl terminal branch (donor at C1 only)."""
    return f"OC[C@H]1O[C@H]{_lab(c1)}[C@H](O)[C@@H](O)[C@H]1O"


def _xylose(c1=None, o4=None, o2=None, o3=None):
    """beta-D-xylopyranose backbone fragment (pentose; no C6)."""
    a = _lab(c1) if c1 is not None else "O"
    b = f"O{_lab(o2)}" if o2 is not None else "O"
    c = f"O{_lab(o3)}" if o3 is not None else "O"
    d = f"O{_lab(o4)}" if o4 is not None else "O"
    return f"C1[C@H]([C@@H]([C@H]([C@@H](O1){a}){b}){c}){d}"


def _araf_donor(c1):
    """alpha-L-arabinofuranosyl terminal branch (donor at C1 only)."""
    return f"C([C@H]1[C@@H]([C@H]([C@@H](O1){_lab(c1)})O)O)O"


def _meglca_donor(c1):
    """4-O-methyl-alpha-D-glucuronic acid terminal branch (donor at C1 only)."""
    return f"CO[C@H]1[C@@H]([C@H]([C@H](O[C@@H]1C(=O)O){_lab(c1)})O)O"


# --------------------------------------------------------------------------------------
def _assert_templates():
    """Gate: every all-free fragment must canonicalize to its verified monomer SMILES."""
    canon = lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s))
    # backbone fragments compare directly to MONOMER_SMILES (all slots free)
    direct = {"bdglc": _hexose("glc"), "bdman": _hexose("man"), "bdxyl": _xylose()}
    for key, frag in direct.items():
        want = canon(molecules.MONOMER_SMILES[key])
        got = canon(frag)
        assert got == want, f"template {key} mismatch:\n  got  {got}\n  want {want}"
    # donor-only fragments: restore the free anomeric OH in place, then compare.
    # gal puts the label where the monomer had a parenthesized (O); xyl/araf/meglca where it had a bare O.
    pairs = [("adgal", _gal_donor(99).replace("%99", "(O)")),
             ("alaraf", _araf_donor(99).replace("%99", "O")),
             ("meglca", _meglca_donor(99).replace("%99", "O"))]
    for key, frag in pairs:
        want = canon(molecules.MONOMER_SMILES[key])
        got = canon(frag)
        assert got == want, f"template {key} mismatch:\n  got  {got}\n  want {want}"


def _finish(smiles, name):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"{name}: assembled SMILES did not parse:\n{smiles}")
    nfrag = len(Chem.GetMolFrags(mol))
    if nfrag != 1:
        raise ValueError(f"{name}: assembled into {nfrag} disconnected pieces (a glycosidic "
                         f"label is unmatched):\n{smiles}")
    return mol


# --------------------------------------------------------------------------------------
# Cellulose: linear beta-1,4 D-glucose homopolymer, reducing end fixed beta.
# Built by the recursive nested-SMILES cassettes (verified to reproduce known cellohexaose).
# --------------------------------------------------------------------------------------
def cellulose(n):
    if not (2 <= n <= 8):
        raise ValueError("cellulose N must be 2..8 (PDB column / embedding limits)")
    FRAG = "OC[C@H]1O[C@@H](O{inner})[C@H](O)[C@@H](O)[C@@H]1O"
    MID = lambda d, inner: f"[C@H]{d}[C@H](O)[C@@H](O)[C@@H](O{inner})O[C@@H]{d}CO"
    CAP = lambda d: f"[C@H]{d}[C@H](O)[C@@H](O)[C@H](O)O[C@@H]{d}CO"   # beta reducing end
    inner = CAP(n)
    for d in range(n - 1, 1, -1):
        inner = MID(d, inner)
    return _finish(FRAG.format(inner=inner), f"cellulose:{n}")


# --------------------------------------------------------------------------------------
# Generic backbone assembler for hemicellulose.
#   units    : list of fragment emitters already bound to their backbone labels
#   We instead build by index here for clarity.
# --------------------------------------------------------------------------------------
def _backbone(seq, branches):
    """
    seq      : list of dicts, one per backbone unit, in NON-reducing -> reducing order:
               {'kind': 'glc'|'man'|'xyl', 'ac2':bool, 'ac3':bool}
    branches : list of (unit_index, position, donor) where position in {'o6','o2','o3'}
               and donor in {'gal','araf','meglca'}.
    Backbone link is unit_i.C1 -(beta 1->4)-> unit_{i+1}.O4.
    Returns the full '.'-joined SMILES.
    """
    n = len(seq)
    label = _L0
    # backbone bond labels: bond k joins unit k (donor C1) to unit k+1 (acceptor O4)
    bond = []
    for _ in range(n - 1):
        bond.append(label); label += 1
    # branch labels + per-unit branch slots
    unit_branch = [dict() for _ in range(n)]
    branch_frags = []
    for (ui, pos, donor) in branches:
        bl = label; label += 1
        unit_branch[ui][pos] = bl
        if donor == "gal":
            branch_frags.append(_gal_donor(bl))
        elif donor == "araf":
            branch_frags.append(_araf_donor(bl))
        elif donor == "meglca":
            branch_frags.append(_meglca_donor(bl))
        else:
            raise ValueError(f"unknown branch donor {donor}")

    frags = []
    for i, u in enumerate(seq):
        c1 = bond[i] if i < n - 1 else None          # donor to next unit
        o4 = bond[i - 1] if i > 0 else None           # acceptor from previous unit
        b = unit_branch[i]
        if u["kind"] in ("glc", "man"):
            frags.append(_hexose(u["kind"], c1=c1, o4=o4,
                                 o6=b.get("o6"), o2=b.get("o2"), o3=b.get("o3"),
                                 ac2=u.get("ac2", False), ac3=u.get("ac3", False)))
        elif u["kind"] == "xyl":
            frags.append(_xylose(c1=c1, o4=o4, o2=b.get("o2"), o3=b.get("o3")))
        else:
            raise ValueError(f"unknown backbone sugar {u['kind']}")
    return ".".join(frags + branch_frags)


def ggm(variant):
    """Softwood galactoglucomannan: beta-1,4 Man/Glc backbone, alpha-1,6 Gal branches."""
    if variant == "compact":
        seq = [{"kind": k} for k in ("man", "man", "glc", "man", "man", "man")]
        branches = [(3, "o6", "gal")]
    elif variant == "extended":
        seq = [{"kind": k} for k in
               ("man", "man", "glc", "man", "man", "man", "glc", "man", "man", "man")]
        seq[1]["ac3"] = True          # native O-acetylation on mannose
        seq[5]["ac2"] = True
        branches = [(3, "o6", "gal"), (7, "o6", "gal")]
    else:
        raise ValueError(variant)
    return _finish(_backbone(seq, branches), f"ggm:{variant}")


def agx(variant):
    """Softwood arabinoglucuronoxylan: beta-1,4 xylan + alpha-1,2 MeGlcA + alpha-1,3 Araf."""
    if variant == "compact":
        seq = [{"kind": "xyl"} for _ in range(5)]
        branches = [(1, "o2", "meglca"), (3, "o3", "araf")]
    elif variant == "extended":
        seq = [{"kind": "xyl"} for _ in range(8)]
        branches = [(1, "o2", "meglca"), (3, "o3", "araf"),
                    (5, "o2", "meglca"), (6, "o3", "araf")]
    else:
        raise ValueError(variant)
    return _finish(_backbone(seq, branches), f"agx:{variant}")


# --------------------------------------------------------------------------------------
# Lignin trimer: a DEFINED representative guaiacyl oligomer (beta-O-4 + beta-5).
# --------------------------------------------------------------------------------------
def lignin_trimer():
    seed = ("COc1cc(C(O)C(CO)Oc2ccc(C3Oc4c(OC)cc(/C=C/CO)cc4C3CO)cc2OC)ccc1O")
    return _finish(seed, "lignin:trimer")


# --------------------------------------------------------------------------------------
_TEMPLATES_OK = False


def build(token):
    """Dispatch a manifest `builder` token to its assembler. Returns an un-embedded Mol."""
    global _TEMPLATES_OK
    if not _TEMPLATES_OK:
        _assert_templates()
        _TEMPLATES_OK = True
    kind, _, arg = token.partition(":")
    if kind == "cellulose":
        return cellulose(int(arg))
    if kind == "ggm":
        return ggm(arg)
    if kind == "agx":
        return agx(arg)
    if kind == "lignin" and arg == "trimer":
        return lignin_trimer()
    raise ValueError(f"unknown builder token {token!r}")


if __name__ == "__main__":
    from rdkit.Chem import rdMolDescriptors
    _assert_templates()
    print("template gate: OK")
    # self-test: cellobiose from the cellulose builder must match known beta-cellobiose
    known_beta_cellobiose = "OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O"
    got = Chem.MolToSmiles(cellulose(2))
    want = Chem.MolToSmiles(Chem.MolFromSmiles(known_beta_cellobiose))
    print("cellobiose self-test:", "OK" if got == want else f"FAIL\n  got {got}\n  want {want}")
    for tok in ("cellulose:6", "ggm:compact", "ggm:extended", "agx:compact",
                "agx:extended", "lignin:trimer"):
        m = build(tok)
        print(f"  {tok:16s} {rdMolDescriptors.CalcMolFormula(m):16s} "
              f"{m.GetNumAtoms()} heavy atoms, frags={len(Chem.GetMolFrags(m))}")
