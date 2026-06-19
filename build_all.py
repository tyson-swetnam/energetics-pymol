"""
build_all.py — one-shot pipeline: manifest -> pdb/ -> scripts/ + renders/ -> contact sheet.

    python build_all.py            # build PDBs, render PNGs, write contact sheet
    python build_all.py --no-png   # build PDBs + .pml only (skip ray-tracing)

Run inside the `energetics-pymol` conda env (see environment.yml / `make env`).
"""
from __future__ import annotations
import os
import sys

import build_pdb
import render
import molecules

ROOT = os.path.dirname(os.path.abspath(__file__))


def contact_sheet():
    """Stitch the per-molecule renders into one catalog PNG grouped by category."""
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError:
        print("(Pillow not installed; skipping contact sheet. `pip install pillow` to enable.)")
        return None

    cols = 4
    cell_w, cell_h, label_h, pad = 320, 240, 26, 8
    title_h = 34
    recs = [r for r in molecules.MANIFEST
            if os.path.exists(os.path.join(ROOT, "renders", r["id"] + ".png"))]

    # group by category, lay out rows
    rows_per_cat = {}
    for cat in molecules.CATEGORIES:
        items = [r for r in recs if r["category"] == cat]
        if items:
            rows_per_cat[cat] = items

    total_rows = sum((len(v) + cols - 1) // cols for v in rows_per_cat.values())
    W = cols * (cell_w + pad) + pad
    H = sum(title_h + ((len(v) + cols - 1) // cols) * (cell_h + label_h + pad) + pad
            for v in rows_per_cat.values()) + pad

    sheet = Image.new("RGB", (W, H), "white")
    draw = ImageDraw.Draw(sheet)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 14)
        tfont = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial Bold.ttf", 20)
    except Exception:
        font = ImageFont.load_default()
        tfont = font

    y = pad
    for cat, items in rows_per_cat.items():
        draw.text((pad, y + 6), cat, fill="black", font=tfont)
        y += title_h
        for i, r in enumerate(items):
            col = i % cols
            if col == 0 and i > 0:
                y += cell_h + label_h + pad
            x = pad + col * (cell_w + pad)
            im = Image.open(os.path.join(ROOT, "renders", r["id"] + ".png")).convert("RGB")
            im.thumbnail((cell_w, cell_h))
            ox = x + (cell_w - im.width) // 2
            sheet.paste(im, (ox, y))
            cap = r["name"]
            if len(cap) > 46:
                cap = cap[:44] + "…"
            draw.text((x, y + cell_h + 4), cap, fill="black", font=font)
        y += cell_h + label_h + pad

    out = os.path.join(ROOT, "renders", "_contact_sheet.png")
    sheet.save(out)
    print(f"contact sheet -> renders/_contact_sheet.png ({W}x{H})")
    return out


def main(argv):
    no_png = "--no-png" in argv
    print("== building PDBs ==")
    rc1 = build_pdb.main(["build_pdb.py"])
    if no_png:
        render.main(["render.py", "--no-png"])
        return rc1
    print("\n== rendering PNGs ==")
    rc2 = render.main(["render.py"])
    print("\n== contact sheet ==")
    contact_sheet()
    return rc1 or rc2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
