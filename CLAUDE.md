# CLAUDE.md

Guidance for working in this repository.

## What this repo is

Reproducible 3D molecular models (PyMOL-ready `.pdb` files + ray-traced figures) of the
organic molecules discussed in the manuscript **"Ecosystems as energy fields"** (Falk,
Swetnam & McKenzie, 2026; `Energetics ms full draft 18 Jun 2026.pdf`). The paper estimates
the embodied chemical-bond energy (enthalpy / higher heating value) stored in forest biomass.

The study system is a **conifer** forest dominated by *Pinus ponderosa* (a softwood /
gymnosperm). **Where chemistry differs by clade, use the softwood representative**:
galactoglucomannan + arabinoglucuronoxylan for hemicellulose, and guaiacyl (G) lignin —
not syringyl-rich hardwood lignin. Cellulose is clade-invariant.

## Environment

A dedicated conda env (`conda-forge`, native Apple Silicon / osx-arm64):

```bash
make env            # conda create -n energetics-pymol -c conda-forge rdkit openbabel pymol-open-source ...
conda activate energetics-pymol
```

- **RDKit is authoritative** for stereochemistry. Open Babel is only a fallback / cross-check;
  `obabel --gen3d` scrambles sugar ring/anomeric stereo, so never use it to (re)build carbohydrates.
- **Never `pip install pymol`** on macOS — only `conda-forge pymol-open-source` has a working
  arm64 build. PyMOL imports as the `pymol` module; rendering runs headless via `pymol -cq`.
- `environment.yml` is the lock (regenerate with `make lock`).

## Repo map

| file | role |
|------|------|
| `molecules.py` | **single source of truth**: the verified molecule manifest (id, name, SMILES/builder, formula, …) |
| `polymers.py`  | oligomer assembly for cellulose / hemicellulose / lignin (returns RDKit Mols) |
| `build_pdb.py` | SMILES/Mol → 3D → validated `pdb/<id>.pdb` |
| `render.py`    | emit `scripts/<id>.pml` + ray-trace `renders/<id>.png` |
| `build_all.py` | orchestrate everything + write `renders/_contact_sheet.png` |
| `build_site.py`| regenerate `docs/molecules.json` + copy assets for the web gallery |
| `Makefile`     | `make env / pdb / render / all / clean` |
| `pdb/ scripts/ renders/` | generated outputs (committed so co-authors don't have to build) |
| `docs/`        | GitHub Pages site (Mol* + 3Dmol gallery); shell is hand-authored, data is generated |

## Web gallery

`docs/` is served by GitHub Pages (`main` → `/docs`) at
https://tyson-swetnam.github.io/energetics-pymol/. The shell (`index.html`, `style.css`,
`app.js`) is hand-authored; run `python build_site.py` to regenerate `docs/molecules.json` and
refresh the `docs/pdb/` + `docs/img/` copies after changing the manifest or rebuilding. The app
serializes Mol* load/clear calls (`state.molstarQueue`) to avoid render-loop races — keep that
if you touch the viewer code. The per-molecule energy gauge values live in `build_site.py`
(`HHV_BY_CATEGORY` / `HHV_OVERRIDE`).

> The unpublished manuscript PDF is **gitignored** — the repo is public; do not commit it.

## How to add or change a molecule

1. Append a record to `MANIFEST` in `molecules.py` (give it an **isomeric** SMILES and a
   3-char `resname`; set `formula` so the build asserts it; set `show_h=True` only for small
   molecules where hydrogens clarify the figure).
2. `make all` (or `python build_pdb.py <id> && python render.py <id>`).
3. Eyeball `renders/<id>.png`.

Use a `builder` token (handled in `polymers.py`) instead of `smiles` only for oligomers.

## Build invariants (do not regress these)

- **Isomeric SMILES everywhere**; `build_pdb.check_stereo` refuses to embed a carbohydrate with
  any undefined ring stereocenter.
- **Deterministic geometry**: ETKDGv3 with a fixed `randomSeed` (+ multi-seed / `useRandomCoords`
  retry). Re-running reproduces identical coordinates.
- **CONECT records are mandatory** (`MolToPDBFile(flavor=0)`). PyMOL needs them for non-standard
  residues; `set connect_mode, 1` in the `.pml` makes PyMOL trust them. Never use `flavor=2`.
- **Oligomers are stereo-safe by construction** (`polymers.py`): glycosidic bonds are SMILES
  ring-closure labels shared across `.`-joined monomer fragments, so no chiral atom's neighbor
  order is ever edited. The template gate + cellobiose self-test (`python polymers.py`) guard this.
- **Lignin / hemicellulose oligomers are *defined representatives***, not the real randomized
  macromolecule — caption figures accordingly. Real lignin is an irregular crosslinked network.
- **Validation gates** (`build_pdb.validate`): molecular formula matches the manifest, CONECT
  present, and the PDB round-trips to the same heavy-atom graph.

## Verify

```bash
python polymers.py          # template gate + cellobiose reconstruction self-test
python build_pdb.py         # all molecules; every line must read [OK]
python render.py --thumb    # fast renders; eyeball for floating fragments / collapsed rings
```
