"""
build_site.py — generate the data + assets for the docs/ GitHub Pages site.

Writes:
    docs/molecules.json   metadata for every molecule (name, formula, SMILES, note, energy)
    docs/pdb/<id>.pdb      copies of the coordinate files (so Pages from /docs can fetch them)
    docs/img/<id>.png      copies of the renders (sidebar thumbnails + viewer poster)
    docs/.nojekyll         serve files verbatim (no Jekyll)

The hand-authored site shell (index.html, style.css, app.js) lives in docs/ and is committed;
this script only (re)generates the data + asset copies so they stay in sync with molecules.py.

    python build_site.py

ENERGY READOUT (the site's signature): each molecule shows an approximate higher heating
value (HHV, MJ/kg) — representative literature values, NOT precise measurements — tying each
structure to the paper's thesis that biomass molecules are energy stores. Lignin and especially
terpene extractives are far more energy-dense than the structural carbohydrates.
"""
from __future__ import annotations
import json
import os
import shutil

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import molecules
import polymers

ROOT = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.join(ROOT, "docs")

# Representative higher heating values (MJ/kg). Sources: paper text (cellulose/hemicellulose
# 15-16; lignin/extractives 25-32), McKendry 2002, and standard fuel/lipid/terpene values.
# These are indicative, for the energy gauge — labelled "approx." in the UI.
HHV_BY_CATEGORY = {
    "Sugar monomers": 15.5,
    "Cellulose": 17.3,
    "Hemicellulose": 16.0,
    "Lignin": 26.5,
    "Fats": 39.5,
    "Extractives": 40.0,
}
HHV_OVERRIDE = {
    "bdglc": 15.6, "glc_open": 15.6, "bdxyl": 15.0, "alaraf": 15.0, "meglca": 14.0,
    "oleic_acid": 39.3, "palmitic_acid": 39.0, "triolein": 39.6, "tripalmitin": 39.1,
    "alpha_pinene": 44.0, "beta_pinene": 44.0, "abietic_acid": 36.5,
}
# Combustion molecules are not fuels:
HHV_KIND = {"o2": "oxidant", "co2": "product", "h2o": "product"}

WHOLE_TREE_HHV = 20.25   # Pinus mean from the paper (Table EC) -> gauge reference line
GAUGE_MAX = 46.0


def _smiles_and_formula(rec):
    if rec["smiles"]:
        mol = Chem.MolFromSmiles(rec["smiles"])
        return rec["smiles"], (rec["formula"] or rdMolDescriptors.CalcMolFormula(mol))
    mol = polymers.build(rec["builder"])
    return Chem.MolToSmiles(mol), (rec["formula"] or rdMolDescriptors.CalcMolFormula(mol))


def build_json():
    mols = []
    for rec in molecules.MANIFEST:
        smiles, formula = _smiles_and_formula(rec)
        kind = HHV_KIND.get(rec["id"], "fuel")
        hhv = None
        if kind == "fuel":
            hhv = HHV_OVERRIDE.get(rec["id"], HHV_BY_CATEGORY.get(rec["category"]))
        mols.append({
            "id": rec["id"],
            "name": rec["name"],
            "resname": rec["resname"],
            "category": rec["category"],
            "formula": formula,
            "smiles": smiles,
            "note": rec["note"],
            "size": rec["size"],
            "showH": rec["show_h"],
            "hhv": hhv,
            "hhvKind": kind,
            "pdb": f"pdb/{rec['id']}.pdb",
            "img": f"img/{rec['id']}.png",
        })
    data = {
        "meta": {
            "title": "Energetics",
            "subtitle": "The molecules that store a forest's energy",
            "paper": "Falk, Swetnam & McKenzie — Ecosystems as energy fields (2026)",
            "wholeTreeHHV": WHOLE_TREE_HHV,
            "gaugeMax": GAUGE_MAX,
            "categories": molecules.CATEGORIES,
            "count": len(mols),
        },
        "molecules": mols,
    }
    os.makedirs(DOCS, exist_ok=True)
    with open(os.path.join(DOCS, "molecules.json"), "w") as fh:
        json.dump(data, fh, indent=2)
    return data


def copy_assets():
    pdb_src, pdb_dst = os.path.join(ROOT, "pdb"), os.path.join(DOCS, "pdb")
    img_src, img_dst = os.path.join(ROOT, "renders"), os.path.join(DOCS, "img")
    os.makedirs(pdb_dst, exist_ok=True)
    os.makedirs(img_dst, exist_ok=True)
    npdb = nimg = 0
    for rec in molecules.MANIFEST:
        p = os.path.join(pdb_src, rec["id"] + ".pdb")
        i = os.path.join(img_src, rec["id"] + ".png")
        if os.path.exists(p):
            shutil.copy2(p, os.path.join(pdb_dst, rec["id"] + ".pdb")); npdb += 1
        if os.path.exists(i):
            shutil.copy2(i, os.path.join(img_dst, rec["id"] + ".png")); nimg += 1
    # serve verbatim (no Jekyll)
    open(os.path.join(DOCS, ".nojekyll"), "w").close()
    return npdb, nimg


def main():
    data = build_json()
    npdb, nimg = copy_assets()
    print(f"docs/molecules.json : {data['meta']['count']} molecules")
    print(f"docs/pdb/           : {npdb} files")
    print(f"docs/img/           : {nimg} files")
    print("docs/.nojekyll      : written")


if __name__ == "__main__":
    main()
