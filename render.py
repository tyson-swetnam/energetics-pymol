"""
render.py — emit a PyMOL .pml script per molecule and ray-trace a publication PNG.

For each manifest record it writes scripts/<id>.pml and (unless --no-png) runs headless
PyMOL (`pymol -cq scripts/<id>.pml`) to ray-trace renders/<id>.png. The .pml files are
standalone: open one in the PyMOL GUI (`@scripts/cellohexaose.pml`) to tweak a figure.

    python render.py                 # all molecules: write .pml + render .png
    python render.py cellohexaose    # one molecule
    python render.py --no-png        # only (re)write the .pml scripts, skip ray-tracing
    python render.py --thumb         # fast low-res thumbnails (verification pass)
"""
from __future__ import annotations
import os
import sys
import subprocess

import molecules

ROOT = os.path.dirname(os.path.abspath(__file__))
PDB_DIR = os.path.join(ROOT, "pdb")
PML_DIR = os.path.join(ROOT, "scripts")
PNG_DIR = os.path.join(ROOT, "renders")

PML_TEMPLATE = """\
# {name}
# {note}
reinitialize
bg_color white
load {pdb}, mol
set connect_mode, 1          # trust CONECT records (non-standard residues)

hide everything
show sticks, mol
{hydro}
set stick_radius, 0.14
set stick_h_scale, 1.0

# grey carbons, heteroatoms by CPK element colour
color grey70, mol and elem C
util.cnc("mol")

set valence, 1               # draw double/aromatic bond lines (lignin rings, C=O, C=C)
set valence_mode, 0

# journal-friendly ray tracing
set ray_opaque_background, 0
set ray_trace_mode, 1        # black outline (cartoon-style edges)
set ray_trace_color, black
set ray_shadows, 0
set antialias, 2
set depth_cue, 0

orient mol
zoom mol, {pad}
ray {w}, {h}
png {png}, dpi={dpi}, ray=1
"""


def pml_for(rec, thumb=False):
    hydro = "" if rec["show_h"] else "hide everything, (mol and hydro)"
    if thumb:
        w, h, dpi, pad = 600, 450, 100, 2.5
    else:
        w, h, dpi, pad = 1600, 1200, 300, 2.0
    return PML_TEMPLATE.format(
        name=rec["name"], note=rec["note"],
        pdb=os.path.join(PDB_DIR, rec["id"] + ".pdb"),
        png=os.path.join(PNG_DIR, rec["id"] + ".png"),
        hydro=hydro, w=w, h=h, dpi=dpi, pad=pad,
    )


def write_pml(rec, thumb=False):
    os.makedirs(PML_DIR, exist_ok=True)
    path = os.path.join(PML_DIR, rec["id"] + ".pml")
    with open(path, "w") as fh:
        fh.write(pml_for(rec, thumb=thumb))
    return path


def render(rec, thumb=False, verbose=True):
    os.makedirs(PNG_DIR, exist_ok=True)
    pml = write_pml(rec, thumb=thumb)
    png = os.path.join(PNG_DIR, rec["id"] + ".png")
    if os.path.exists(png):
        os.remove(png)
    proc = subprocess.run(["pymol", "-cq", pml],
                          capture_output=True, text=True, cwd=ROOT)
    ok = os.path.exists(png) and os.path.getsize(png) > 0
    if verbose:
        status = "OK" if ok else "FAIL"
        print(f"[{status}] {rec['id']:16s} -> renders/{rec['id']}.png")
        if not ok and proc.stderr:
            print("       " + proc.stderr.strip().splitlines()[-1])
    return ok


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    flags = {a for a in argv[1:] if a.startswith("--")}
    no_png = "--no-png" in flags
    thumb = "--thumb" in flags

    recs = [molecules.BY_ID[i] for i in args] if args else molecules.MANIFEST
    if no_png:
        for r in recs:
            write_pml(r, thumb=thumb)
        print(f"Wrote {len(recs)} .pml script(s) to scripts/ (no rendering).")
        return 0

    results = [render(r, thumb=thumb) for r in recs]
    nfail = sum(1 for ok in results if not ok)
    print(f"\nRendered {len(results)} molecule(s); {nfail} failed.")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
