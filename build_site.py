"""
build_site.py — generate the data + assets for the docs/ GitHub Pages site.

Writes:
    docs/molecules.json   metadata for every molecule (name, formula, SMILES, note, energy)
    docs/pdb/<id>.pdb      copies of the coordinate files (so Pages from /docs can fetch them)
    docs/img/<id>.png      copies of the renders (sidebar thumbnails + viewer poster)
    docs/.nojekyll         serve files verbatim (no Jekyll)

and, for search engines and LLM agents (all regenerated, all deterministic — no timestamps):
    docs/robots.txt        crawl policy + pointers to the sitemap and llms.txt
    docs/sitemap.xml       every page + per-molecule PDB/PNG
    docs/llms.txt          llmstxt.org index: what the site is, where the machine-readable data is
    docs/llms-full.txt     the whole site as one Markdown file (catalog + ontology + learn-page text)
    docs/molecules.jsonld  Bioschemas/schema.org Dataset of MolecularEntity records (ChEBI sameAs)

The hand-authored site shell (index.html, style.css, app.js) lives in docs/ and is committed;
this script only (re)generates the data + asset copies so they stay in sync with molecules.py.

    python build_site.py

ENERGY READOUT (the site's signature): each molecule shows an approximate higher heating
value (HHV, MJ/kg) — representative literature values, NOT precise measurements — tying each
structure to the paper's thesis that biomass molecules are energy stores. Lignin and especially
terpene extractives are far more energy-dense than the structural carbohydrates.
"""
from __future__ import annotations
import html
import html.parser
import json
import os
import re
import shutil

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import molecules
import ontology
import polymers

ROOT = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.join(ROOT, "docs")
SITE = "https://tyson-swetnam.github.io/energetics-pymol/"
REPO = "https://github.com/tyson-swetnam/energetics-pymol"
MESA_PROJECT = "/iplant/home/tswetnam/energetics-pymol"   # CyVerse Data Store (iRODS) + DuckLake

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
    "Fossil carbon": 32.0,   # energy-dense reduced carbon (dry/model basis); approx.
    "Intermediates": 22.0,   # activated carriers — shown for reference, not energy stores
}
HHV_OVERRIDE = {
    "bdglc": 15.6, "glc_open": 15.6, "bdxyl": 15.0, "alaraf": 15.0, "meglca": 14.0,
    "oleic_acid": 39.3, "palmitic_acid": 39.0, "triolein": 39.6, "tripalmitin": 39.1,
    "alpha_pinene": 44.0, "beta_pinene": 44.0, "abietic_acid": 36.5,
    # fossil carbon (approx., dry/dmmf or model basis; coalification concentrates energy)
    "charcoal": 33.0, "lignite": 24.0, "sub_bituminous": 28.0, "bituminous": 33.0,
    "anthracite": 34.0, "kerogen": 30.0, "crude_oil": 44.0,
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


def _with_iri(terms):
    return [{**t, "iri": ontology.iri(t["curie"])} for t in terms]


def build_json():
    mols = []
    for rec in molecules.MANIFEST:
        smiles, formula = _smiles_and_formula(rec)
        inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
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
            "inchikey": inchikey,
            "note": rec["note"],
            "size": rec["size"],
            "showH": rec["show_h"],
            "hhv": hhv,
            "hhvKind": kind,
            "pdb": f"pdb/{rec['id']}.pdb",
            "img": f"img/{rec['id']}.png",
            "ontology": {k: v for k, v in ontology.terms_for(rec["id"]).items()
                         if k in ("structure", "represents")},   # shared layers live in meta
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
            "site": SITE,
            "repo": REPO,
            "mesaProject": MESA_PROJECT,
            "ontology": {
                "matchLevels": ontology.MATCH_LEVELS,
                "project": _with_iri(ontology.PROJECT_TERMS),
                "categories": {c: _with_iri(ts) for c, ts in ontology.CATEGORY_TERMS.items()},
            },
        },
        "molecules": mols,
    }
    os.makedirs(DOCS, exist_ok=True)
    with open(os.path.join(DOCS, "molecules.json"), "w", encoding="utf-8", newline="\n") as fh:
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


# --------------------------------------------------------------------------------------
# Agent-readable files: robots.txt, sitemap.xml, llms.txt, llms-full.txt, molecules.jsonld
# --------------------------------------------------------------------------------------
PAGES = [
    ("index.html", "Molecule gallery",
     "Interactive 3D viewer (Mol* / 3Dmol.js) for every molecule, with an approximate "
     "higher-heating-value energy gauge."),
    ("learn.html", "How plants build it — and how nature releases it",
     "Narrative explainer: photosynthesis, the biosynthesis of cellulose, hemicellulose, lignin, "
     "fats and resin, then fire, consumption, decay and charcoal/fossil carbon."),
]

# Crawlers that identify as AI/LLM agents; listed explicitly so the allow is unambiguous.
AI_AGENTS = ["GPTBot", "ChatGPT-User", "OAI-SearchBot", "ClaudeBot", "Claude-User",
             "Claude-SearchBot", "anthropic-ai", "PerplexityBot", "Perplexity-User",
             "Google-Extended", "Applebot-Extended", "CCBot", "Meta-ExternalAgent"]


class _TextExtractor(html.parser.HTMLParser):
    """Minimal HTML -> Markdown-ish text for llms-full.txt (headings, paragraphs, lists).

    Reaction-chain lists (`ol.rxn`, `ol.pathway`) collapse to one line, e.g.
    `CO₂ (from air) →[RuBisCO + RuBP, fixation]→ 3-phosphoglycerate ...`.
    """
    SKIP = {"script", "style", "svg", "noscript", "header", "nav", "button"}
    BLOCK = {"p", "div", "section", "figure", "figcaption", "details", "summary", "ol", "ul",
             "br", "table", "tr", "blockquote", "footer"}
    HEAD = {"h1": "#", "h2": "##", "h3": "###", "h4": "####"}
    INLINE_LISTS = {"rxn", "pathway"}
    SEGMENTS = {"seg-cell": "cellulose", "seg-hemi": "hemicellulose", "seg-lig": "lignin",
                "seg-ext": "extractives", "seg-oth": "other"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.out, self.skip, self.inline, self.li = [], 0, 0, []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        cls = set((a.get("class") or "").split())
        style = dict(kv.split(":", 1) for kv in (a.get("style") or "").replace(" ", "").split(";") if ":" in kv)
        if tag in self.SKIP:
            self.skip += 1
        elif self.skip:
            return
        elif tag in ("ol", "ul") and cls & self.INLINE_LISTS:
            self.inline += 1
            self.out.append("\n\n")
        elif self.inline and tag == "li":
            kind = "step" if "step" in cls else "head" if "head" in cls else "met"
            self.li.append(kind)
            self.out.append({"step": " →[", "head": "\n\n**", "met": " "}[kind])
        elif tag == "small":
            self.out.append(" (")
        elif tag == "span" and "cof" in cls:     # enzyme | cofactor pair inside a pathway step
            self.out.append("; ")
        elif tag == "span" and "comp-lab" in cls:   # one composition row per tree part
            self.out.append("\n- ")
        elif tag == "i" and cls & set(self.SEGMENTS) and "width" in style:
            # composition bar: the percentages exist only as CSS widths
            self.out.append(f" {self.SEGMENTS[(cls & set(self.SEGMENTS)).pop()]} {style['width']},")
        elif tag == "span" and "thermo-bar" in cls:
            # pyrolysis bar: CSS left/width on a 100-900 °C axis -> temperature range
            lo, w = float(style["left"].rstrip("%")), float(style["width"].rstrip("%"))
            t = lambda pct: int(round((100 + 8 * pct) / 10.0) * 10)
            self.out.append(f" (~{t(lo)}–{t(lo + w)} °C)")
        elif tag in ("span", "b") and self.out and not self.out[-1][-1:].isspace():
            self.out.append(" ")                 # adjacent label/value cells, axis ticks
        elif tag in self.HEAD:
            self.out.append(f"\n\n{self.HEAD[tag]} ")
        elif tag == "li":
            self.out.append("\n- ")
        elif tag == "summary":
            self.out.append("\n\n**")
        elif tag in self.BLOCK:
            self.out.append("\n\n")

    def handle_endtag(self, tag):
        if tag in self.SKIP:
            self.skip = max(0, self.skip - 1)
        elif self.skip:
            return
        elif self.inline and tag in ("ol", "ul"):
            self.inline -= 1
            self.out.append("\n\n")
        elif self.inline and tag == "li":
            kind = self.li.pop() if self.li else "met"
            self.out.append({"step": "]→", "head": ":**", "met": ""}[kind])
        elif tag == "small":
            self.out.append(")")
        elif tag == "summary":
            self.out.append("**\n")
        elif tag in self.BLOCK or tag in self.HEAD:
            self.out.append("\n\n")

    def handle_data(self, data):
        if not self.skip:
            self.out.append(re.sub(r"\s+", " ", data))

    def text(self):
        t = "".join(self.out)
        t = re.sub(r",(\s*\n)", r"\1", t)          # trailing comma after the last bar segment
        t = re.sub(r"\( +", "(", t)
        t = re.sub(r" +\)", ")", t)
        t = re.sub(r"\[ +", "[", t)
        t = re.sub(r" +\]", "]", t)
        t = re.sub(r" {2,}", " ", t)
        t = re.sub(r"\n- \s*\n", "\n", t)
        lines = [ln.rstrip() if ln.startswith("- ") else ln.strip() for ln in t.splitlines()]
        return re.sub(r"\n{3,}", "\n\n", "\n".join(lines)).strip()


def _page_text(name):
    with open(os.path.join(DOCS, name), encoding="utf-8") as fh:
        src = fh.read()
    main_html = re.search(r"<main\b.*?</main>", src, re.S)
    p = _TextExtractor()
    p.feed(main_html.group(0) if main_html else src)
    return p.text()


def _term_md(t):
    rel = f" — {t['relation']}" if t.get("relation") else ""
    return f"[{t['curie']}]({t['iri']}) {t['label']}{rel}"


def write_robots():
    lines = ["# robots.txt — Energetics (energetics-pymol), a GitHub Pages project site.",
             "# Everything here is public, static and safe to crawl, index and use for AI",
             "# retrieval. Machine-readable entry points for LLM agents:",
             f"#   {SITE}llms.txt       (index)",
             f"#   {SITE}llms-full.txt  (entire site as one Markdown file)",
             f"#   {SITE}molecules.json (structured catalog incl. OBO ontology terms)",
             "",
             "User-agent: *", "Allow: /", ""]
    for a in AI_AGENTS:
        lines += [f"User-agent: {a}"]
    lines += ["Allow: /", "", f"Sitemap: {SITE}sitemap.xml", ""]
    with open(os.path.join(DOCS, "robots.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines))


def write_sitemap(data):
    urls = [SITE] + [SITE + p for p, _, _ in PAGES if p != "index.html"]
    urls += [SITE + f for f in ("llms.txt", "llms-full.txt", "molecules.json", "molecules.jsonld")]
    for m in data["molecules"]:
        urls += [SITE + m["pdb"], SITE + m["img"]]
    body = "\n".join(f"  <url><loc>{html.escape(u)}</loc></url>" for u in urls)
    with open(os.path.join(DOCS, "sitemap.xml"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write('<?xml version="1.0" encoding="UTF-8"?>\n'
                 '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
                 f"{body}\n</urlset>\n")


def _intro(data):
    meta = data["meta"]
    return [
        f"> {meta['count']} reproducible 3D molecular models (PDB with CONECT bonds + ray-traced "
        "PNG) of the organic molecules that store chemical-bond energy in a Pinus ponderosa "
        f"(ponderosa pine, conifer/softwood) forest, built for the manuscript \"{meta['paper']}\". "
        "Every model is annotated with OBO Foundry ontology terms (ChEBI, GO, ENVO, PO, NCBITaxon).",
        "",
        "Key facts for agents:",
        "",
        "- Study system is a conifer, so hemicellulose = softwood galactoglucomannan + "
        "arabinoglucuronoxylan and lignin = guaiacyl (G) lignin. Cellulose is clade-invariant.",
        "- Oligomers (cellulose N=2–8, hemicellulose fragments, lignin dimers/trimer) and the "
        "fossil-carbon entries are *defined representative models*, not the real randomized "
        "macromolecules or materials; cite them as such.",
        "- Energy values (HHV, MJ/kg) are approximate literature values for the gauge, not "
        f"measurements. Whole-tree Pinus reference: {meta['wholeTreeHHV']} MJ/kg.",
        "- Ontology `match` levels: " + "; ".join(f"`{k}` = {v}" for k, v in
                                                   ontology.MATCH_LEVELS.items()) + ".",
        f"- Source code and build pipeline: {REPO}. Data mirror (CyVerse Data Store, "
        f"MESA/DuckLake-tracked metadata): iRODS `{MESA_PROJECT}`.",
    ]


def write_llms_txt(data):
    out = [f"# {data['meta']['title']} — {data['meta']['subtitle']}", ""] + _intro(data) + [""]
    out += ["## Pages", ""]
    out += [f"- [{t}]({SITE}{p}): {d}" for p, t, d in PAGES]
    out += ["", "## Machine-readable data", "",
            f"- [llms-full.txt]({SITE}llms-full.txt): the entire site — full molecule catalog "
            "with SMILES, InChIKeys, formulas, energy values and ontology terms, plus the full "
            "text of the explainer page — as one Markdown file",
            f"- [molecules.json]({SITE}molecules.json): structured catalog (JSON), one record per "
            "molecule with `ontology.structure` / `ontology.represents`; category- and "
            "project-level terms are in `meta.ontology`",
            f"- [molecules.jsonld]({SITE}molecules.jsonld): schema.org Dataset of Bioschemas "
            "MolecularEntity records cross-referenced to ChEBI (JSON-LD; `sameAs` only for exact matches)",
            f"- [sitemap.xml]({SITE}sitemap.xml): every page and file", "",
            "## Molecules (3D coordinate files)", ""]
    for cat in data["meta"]["categories"]:
        out += [f"### {cat}", ""]
        for m in (x for x in data["molecules"] if x["category"] == cat):
            s = m["ontology"]["structure"]
            out.append(f"- [{m['name']}]({SITE}{m['pdb']}): {m['formula']}; "
                       f"{s['curie']} ({s['match']})")
        out.append("")
    out += ["## Optional", "",
            f"- [Source repository]({REPO}): RDKit/PyMOL build pipeline (`molecules.py` is the "
            "single source of truth; `ontology.py` holds the OBO term tables)",
            f"- [Rendered figures]({SITE}img/cellohexaose.png): one ray-traced PNG per molecule at "
            "`img/<id>.png` (each molecules.json record's `img` field)", ""]
    with open(os.path.join(DOCS, "llms.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out))


def write_llms_full(data):
    meta = data["meta"]
    out = [f"# {meta['title']} — {meta['subtitle']}", "",
           f"Source: {SITE} · Repository: {REPO}", ""] + _intro(data) + [""]
    onto = meta["ontology"]
    proj = onto["project"]
    out += ["## Project-level ontology terms", ""] + [f"- {_term_md(t)}" for t in proj] + [""]
    out += ["## Molecule catalog", ""]
    for cat in meta["categories"]:
        mols = [x for x in data["molecules"] if x["category"] == cat]
        out += [f"### {cat}", ""]
        if mols:
            cterms = onto["categories"].get(cat, [])
            out += ["Category terms: " + ("; ".join(_term_md(t) for t in cterms) or
                                          "none shared — see each molecule"), ""]
        for m in mols:
            o = m["ontology"]
            energy = (f"~{m['hhv']} MJ/kg (approx. HHV)" if m["hhv"] is not None
                      else f"n/a ({m['hhvKind']} — not a fuel)")
            out += [f"#### {m['name']}", "",
                    f"- id: `{m['id']}` · residue `{m['resname']}` · formula {m['formula']}",
                    f"- SMILES (isomeric): `{m['smiles']}`",
                    f"- InChIKey: `{m['inchikey']}`",
                    f"- energy: {energy}",
                    f"- ChEBI structure: {_term_md(o['structure'])} — match `{o['structure']['match']}`"]
            out += [f"- represents: {_term_md(t)}" for t in o["represents"]]
            out += [f"- note: {m['note']}",
                    f"- files: [PDB]({SITE}{m['pdb']}) · [PNG]({SITE}{m['img']}) · "
                    f"[3D viewer]({SITE}#{m['id']})", ""]
    for page, title, _ in PAGES:
        if page == "index.html":
            continue
        out += ["---", "", f"## Page: {title}", "", f"Source: {SITE}{page}", "", _page_text(page), ""]
    with open(os.path.join(DOCS, "llms-full.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(out))


IS_A = ("is a", "is a member of the class")   # relations that are true rdf:type claims


def _dataset_node(data):
    meta = data["meta"]
    return {
        "@type": "Dataset",
        "@id": SITE,
        "dct:conformsTo": {"@id": "https://bioschemas.org/profiles/Dataset/1.0-RELEASE"},
        "name": f"{meta['title']} — {meta['subtitle']}",
        "description": _intro(data)[0].lstrip("> "),
        "url": SITE,
        "sameAs": {"@id": REPO},
        "citation": meta["paper"],
        "keywords": [t["label"] for t in meta["ontology"]["project"]] + meta["categories"],
        "about": [{"@id": t["iri"], "name": t["label"]} for t in meta["ontology"]["project"]],
        "distribution": [
            {"@type": "DataDownload", "contentUrl": SITE + "molecules.jsonld",
             "encodingFormat": "application/ld+json"},
            {"@type": "DataDownload", "contentUrl": SITE + "molecules.json",
             "encodingFormat": "application/json"},
            {"@type": "DataDownload", "contentUrl": SITE + "llms-full.txt",
             "encodingFormat": "text/markdown"}],
    }


JSONLD_CONTEXT = ["https://schema.org/", {"dct": "http://purl.org/dc/terms/"}]


def write_jsonld(data):
    """schema.org/Bioschemas graph: one Dataset node + one MolecularEntity node per molecule.
    Only exact InChIKey matches become `sameAs`; looser ChEBI matches are a qualified
    `identifier`. Only true is-a relations become `additionalType`; ChEBI roles become
    `chemicalRole`. The richer relations (precursor of, produced by, ...) stay in molecules.json."""
    dataset = _dataset_node(data)
    ents = []
    for m in data["molecules"]:
        o, s = m["ontology"], m["ontology"]["structure"]
        ent = {
            "@type": "MolecularEntity",
            "@id": f"{SITE}#{m['id']}",
            "dct:conformsTo": {"@id": "https://bioschemas.org/profiles/MolecularEntity/0.5-RELEASE"},
            "identifier": [m["id"], {"@type": "PropertyValue", "propertyID": "ChEBI",
                                     "name": f"ChEBI structure ({s['match']} match)",
                                     "value": s["curie"], "url": s["iri"]}],
            "name": m["name"],
            "description": m["note"],
            "molecularFormula": m["formula"],
            "smiles": m["smiles"],
            "inChIKey": m["inchikey"],
            "url": f"{SITE}#{m['id']}",
            "image": SITE + m["img"],
            "hasRepresentation": SITE + m["pdb"],
            "subjectOf": [{"@id": SITE},
                          {"@type": "DataDownload", "contentUrl": SITE + m["pdb"],
                           "encodingFormat": "chemical/x-pdb"}],
        }
        if s["match"] == "exact":
            ent["sameAs"] = {"@id": s["iri"]}
        types = [{"@id": t["iri"]} for t in o["represents"] if t["relation"] in IS_A]
        roles = [{"@type": "DefinedTerm", "@id": t["iri"], "name": t["label"],
                  "inDefinedTermSet": "ChEBI"} for t in o["represents"] if t["relation"] == "has role"]
        if types:
            ent["additionalType"] = types
        if roles:
            ent["chemicalRole"] = roles
        ents.append(ent)
    dataset["about"] += [{"@id": e["@id"]} for e in ents]
    doc = {"@context": JSONLD_CONTEXT, "@graph": [dataset] + ents}
    with open(os.path.join(DOCS, "molecules.jsonld"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(doc, fh, indent=2, ensure_ascii=False)
        fh.write("\n")


JSONLD_BEGIN, JSONLD_END = "<!-- build_site:dataset-jsonld -->", "<!-- /build_site:dataset-jsonld -->"


def embed_dataset_jsonld(data):
    """Inline the Dataset node in index.html (harvesters read embedded JSON-LD, not <link>s)."""
    path = os.path.join(DOCS, "index.html")
    with open(path, encoding="utf-8") as fh:
        html_src = fh.read()
    node = {"@context": JSONLD_CONTEXT, **_dataset_node(data)}
    block = (f"{JSONLD_BEGIN}\n<script type=\"application/ld+json\">\n"
             f"{json.dumps(node, indent=2, ensure_ascii=False)}\n</script>\n{JSONLD_END}")
    if JSONLD_BEGIN in html_src:
        pre, rest = html_src.split(JSONLD_BEGIN, 1)
        html_src = pre + block + rest.split(JSONLD_END, 1)[1]
    else:
        html_src = html_src.replace("</head>", block + "\n</head>", 1)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(html_src)


def write_agent_files(data):
    write_robots()
    write_sitemap(data)
    write_llms_txt(data)
    write_llms_full(data)
    write_jsonld(data)
    embed_dataset_jsonld(data)


def main():
    data = build_json()
    npdb, nimg = copy_assets()
    write_agent_files(data)
    print(f"docs/molecules.json : {data['meta']['count']} molecules")
    print(f"docs/pdb/           : {npdb} files")
    print(f"docs/img/           : {nimg} files")
    print("docs/.nojekyll      : written")
    print("docs/robots.txt, sitemap.xml, llms.txt, llms-full.txt, molecules.jsonld : written")


if __name__ == "__main__":
    main()
