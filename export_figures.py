"""
export_figures.py — high-resolution, transparent-background molecule PNGs for figure work
(e.g. recreating the manuscript's Figure A — cellulose, hemicellulose, lignin — in Claude Design).

Renders to export/<id>.png at high resolution (default 3000x2250, 600 dpi metadata), transparent
background, ball-and-stick (grey carbons + CPK heteroatoms, black outline). Drop straight onto any
canvas — the alpha channel means no white box to mask out.

    python export_figures.py                 # the Figure-A set (3 polymers + building blocks)
    python export_figures.py --all           # every molecule in the manifest
    python export_figures.py cellohexaose lignin_trimer   # named ids
    python export_figures.py --size 4000     # override long-edge pixels

Run inside the `energetics-pymol` conda env (needs pymol-open-source).
"""
from __future__ import annotations
import os
import sys
import tempfile
import subprocess

import molecules

ROOT = os.path.dirname(os.path.abspath(__file__))
PDB_DIR = os.path.join(ROOT, "pdb")
EXPORT_DIR = os.path.join(ROOT, "export")

# The molecules that actually appear in the PDF's molecular figure (Figure A), with strong
# representatives + the building blocks you'd use to compose a "cellulose / hemicellulose / lignin"
# figure. Use --all for the full 33-molecule manifest.
FIGURE_SET = [
    "cellohexaose", "cellooctaose",                 # cellulose
    "ggm_compact", "ggm_extended", "agx_compact",   # hemicellulose (softwood)
    "lignin_trimer", "lignin_bo4", "lignin_b5", "lignin_bb",  # lignin
    "bdglc", "coniferyl",                           # building blocks (glucose, monolignol)
]

PML_TEMPLATE = """\
# {name}
reinitialize
bg_color white
load {pdb}, mol
set connect_mode, 1
hide everything
show sticks, mol
show spheres, mol
{hydro}
set stick_radius, 0.15
set sphere_scale, 0.23
set stick_h_scale, 1.0

color grey70, mol and elem C
util.cnc("mol")

set valence, 1
set valence_mode, 0

# transparent, journal-style ray trace
set ray_opaque_background, 0
set ray_trace_mode, 1
set ray_trace_color, black
set ray_shadows, 0
set antialias, 2
set depth_cue, 0

orient mol
zoom mol, 2.5
# size is set ON the png command (ray=1 ray-traces at that exact pixel size);
# a separate `ray W,H` would be overridden by png's own re-trace at viewport size.
png {png}, width={w}, height={h}, dpi={dpi}, ray=1
"""


def pml_for(rec, w, h, dpi):
    hydro = "" if rec["show_h"] else "hide everything, (mol and hydro)"
    return PML_TEMPLATE.format(
        name=rec["name"],
        pdb=os.path.join(PDB_DIR, rec["id"] + ".pdb"),
        png=os.path.join(EXPORT_DIR, rec["id"] + ".png"),
        hydro=hydro, w=w, h=h, dpi=dpi,
    )


def export(rec, w, h, dpi, verbose=True):
    os.makedirs(EXPORT_DIR, exist_ok=True)
    png = os.path.join(EXPORT_DIR, rec["id"] + ".png")
    if os.path.exists(png):
        os.remove(png)
    with tempfile.NamedTemporaryFile("w", suffix=".pml", delete=False) as fh:
        fh.write(pml_for(rec, w, h, dpi))
        pml = fh.name
    try:
        proc = subprocess.run(["pymol", "-cq", pml], capture_output=True, text=True, cwd=ROOT)
    finally:
        os.unlink(pml)
    ok = os.path.exists(png) and os.path.getsize(png) > 0
    if verbose:
        size = f"{os.path.getsize(png)//1024} KB" if ok else "-"
        print(f"[{'OK' if ok else 'FAIL'}] {rec['id']:16s} -> export/{rec['id']}.png  ({size})")
        if not ok and proc.stderr:
            print("       " + proc.stderr.strip().splitlines()[-1])
    return ok


def main(argv):
    flags = {a for a in argv[1:] if a.startswith("--")}
    args = [a for a in argv[1:] if not a.startswith("--")]

    # --size N sets the long edge; height keeps a 4:3 frame
    long_edge = 3000
    for i, a in enumerate(argv):
        if a == "--size" and i + 1 < len(argv):
            long_edge = int(argv[i + 1]); args = [x for x in args if x != argv[i + 1]]
    w, h, dpi = long_edge, int(long_edge * 0.75), 600

    if "--all" in flags:
        recs = molecules.MANIFEST
    elif args:
        recs = [molecules.BY_ID[i] for i in args]
    else:
        recs = [molecules.BY_ID[i] for i in FIGURE_SET]

    print(f"Exporting {len(recs)} molecule(s) at {w}x{h}, transparent, to export/\n")
    results = [export(r, w, h, dpi) for r in recs]
    nfail = sum(1 for ok in results if not ok)
    print(f"\nExported {len(results) - nfail}/{len(results)} to export/  (transparent PNG, {w}x{h}).")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
