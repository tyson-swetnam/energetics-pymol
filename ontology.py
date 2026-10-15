"""
ontology.py — OBO Foundry ontology annotations for every molecule in the manifest.

Each model is tagged with terms from OBO Foundry ontologies, resolved and cross-checked
against the EMBL-EBI Ontology Lookup Service (OLS4) by an adversarial verification pass:

    CHEBI      Chemical Entities of Biological Interest  (what the modeled structure IS)
    GO         Gene Ontology                             (the process that makes it)
    ENVO       Environment Ontology                      (biome / material context)
    PO         Plant Ontology                            (the plant tissue it lives in)
    NCBITaxon  NCBI organismal classification            (Pinus ponderosa)

Three layers, most specific first:

    MOLECULE_TERMS[id]["structure"]   the ChEBI entity for the modeled structure, with an
                                      honest `match` level (see MATCH_LEVELS)
    MOLECULE_TERMS[id]["represents"]  what the model stands for (polymer, material, role)
    CATEGORY_TERMS[category]          shared by every molecule in a manifest category
    PROJECT_TERMS                     shared by the whole collection

The annotations are written into each PDB as REMARK 950 records (an unassigned remark
number in wwPDB v3.3, so no reader confuses them with standard remarks) plus a COMPND
name. RDKit, PyMOL, Mol* and 3Dmol all ignore them, so coordinates/CONECT are untouched.

    python ontology.py              # print the annotation table
    python ontology.py --annotate   # (re)write the REMARK 950 block in pdb/*.pdb in place
    python ontology.py --check      # gate: every manifest id annotated + every PDB carries it
    python ontology.py --verify-ols # network: labels/obsolescence + InChIKey vs model, live OLS4
"""
from __future__ import annotations

import os
import sys

import molecules

ROOT = os.path.dirname(os.path.abspath(__file__))
PDB_DIR = os.path.join(ROOT, "pdb")

OBO = "http://purl.obolibrary.org/obo/"
REMARK = "REMARK 950 "          # 11 cols -> 69 cols of payload per line
PAYLOAD = 80 - len(REMARK)

MATCH_LEVELS = {
    "exact": "same full InChIKey as the ChEBI entity",
    "stereo_unspecified": "same constitution (InChIKey block 1); stereo differs or is undefined",
    "close": "near-identical entity (e.g. anomer/ionization/tautomer variant)",
    "broader": "ChEBI parent class; no entity exists for this exact model",
    "none": "no suitable ChEBI entity",
}


def iri(curie: str) -> str:
    """OBO PURL for a CURIE, e.g. CHEBI:15903 -> http://purl.obolibrary.org/obo/CHEBI_15903."""
    prefix, local = curie.split(":", 1)
    return f"{OBO}{prefix}_{local}"


def _t(curie, label, relation=""):
    return dict(curie=curie, label=label, relation=relation)


# --------------------------------------------------------------------------------------
# Verified term tables (OLS4-confirmed; see module docstring)
# --------------------------------------------------------------------------------------
def _m(curie, label, match, represents):
    return dict(structure=dict(curie=curie, label=label, match=match), represents=represents)


PROJECT_TERMS: list[dict] = [
    _t('NCBITaxon:55062', 'Pinus ponderosa', 'source organism (study system)'),
    _t('ENVO:01000211', 'temperate coniferous forest biome', 'biome'),
    _t('ENVO:01001243', 'forest ecosystem', 'ecosystem whose bond energy is estimated'),
    _t('PO:0005848', 'secondary xylem', 'plant tissue (wood) holding the polymers'),
    _t('ENVO:00002040', 'wood', 'environmental material'),
    _t('ENVO:01000155', 'organic material', 'biomass (closest OBO class)'),
    _t('GO:0015979', 'photosynthesis', 'process that captures the stored energy'),
]

CATEGORY_TERMS: dict[str, list[dict]] = {
    'Sugar monomers': [
        _t('CHEBI:35381', 'monosaccharide', 'chemical class'),
        _t('GO:0046364', 'monosaccharide biosynthetic process', 'produced by'),
    ],
    'Cellulose': [
        _t('CHEBI:18246', '(1->4)-beta-D-glucan', 'cellulose; models are oligomers of it'),
        _t('GO:0030244', 'cellulose biosynthetic process', 'produced by'),
        _t('GO:0052324', 'plant-type cell wall cellulose biosynthetic process', 'produced by'),
    ],
    'Hemicellulose': [
        _t('CHEBI:61266', 'hemicellulose', 'defined representative fragment of'),
        _t('GO:0070592', 'cell wall polysaccharide biosynthetic process', 'produced by'),
    ],
    'Lignin': [
        _t('CHEBI:6457', 'lignin', 'monomer precursor or linkage model of'),
        _t('GO:0009809', 'lignin biosynthetic process', 'produced by / participant in'),
    ],
    'Combustion': [
        _t('ENVO:01000839', 'combustion process', 'reactant or product'),
        _t('ENVO:01000791', 'forest fire', 'reactant or product'),
        _t('GO:0009060', 'aerobic respiration', 'same reactants/products'),
    ],
    'Fats': [
        _t('CHEBI:18059', 'lipid', 'chemical class'),
        _t('GO:0006633', 'fatty acid biosynthetic process', 'acyl chains produced by'),
    ],
    'Extractives': [
        _t('CHEBI:24913', 'isoprenoid', 'chemical class'),
        _t('GO:0016114', 'terpenoid biosynthetic process', 'produced by'),
    ],
    'Fossil carbon': [],
    'Intermediates': [],
}

MOLECULE_TERMS: dict[str, dict] = {
    'bdglc': _m('CHEBI:15903', 'beta-D-glucose', 'exact', [
        _t('CHEBI:18246', '(1->4)-beta-D-glucan', 'is the repeat unit (monomer) of'),
        _t('CHEBI:17634', 'D-glucose', 'is the beta-pyranose anomer of'),
        _t('CHEBI:4167', 'D-glucopyranose', 'is a'),
    ]),
    'glc_open': _m('CHEBI:42758', 'aldehydo-D-glucose', 'exact', [
        _t('CHEBI:17634', 'D-glucose', 'is an open-chain tautomer of'),
        _t('CHEBI:33917', 'aldohexose', 'is a'),
    ]),
    'bdman': _m('CHEBI:28563', 'beta-D-mannose', 'exact', [
        _t('CHEBI:24144', 'galactoglucomannan', 'is the principal backbone monomer unit of'),
        _t('CHEBI:17020', 'glucomannan', 'is a monomer unit of'),
        _t('CHEBI:61266', 'hemicellulose', 'is a constituent sugar of softwood'),
    ]),
    'bdgal': _m('CHEBI:27667', 'beta-D-galactose', 'exact', [
        _t('CHEBI:24144', 'galactoglucomannan', 'is a monomer of'),
        _t('CHEBI:61266', 'hemicellulose', 'is a constituent sugar of'),
    ]),
    'bdxyl': _m('CHEBI:28161', 'beta-D-xylose', 'exact', [
        _t('CHEBI:37166', 'xylan', 'is the monomer unit of'),
        _t('CHEBI:60938', 'glucuronoxylan', 'is the backbone monomer of'),
        _t('CHEBI:61266', 'hemicellulose', 'is a constituent monomer of'),
    ]),
    'alaraf': _m('CHEBI:28772', 'alpha-L-arabinofuranose', 'exact', [
        _t('CHEBI:24300', 'glucuronoarabinoxylans', 'is a side-branch residue of'),
        _t('CHEBI:22603', 'arabinoxylans', 'is a side-branch residue of'),
        _t('CHEBI:61266', 'hemicellulose', 'is a monomer constituent of'),
        _t('CHEBI:30849', 'L-arabinose', 'is the alpha-furanose ring form of'),
    ]),
    'meglca': _m('CHEBI:234000', '4-O-methyl-beta-D-glucuronate', 'close', [
        _t('CHEBI:17856', 'glucuronoxylan 4-O-methyl-D-glucuronate',
           'is the side-chain residue of'),
        _t('CHEBI:60938', 'glucuronoxylan', 'is a substituent monomer of'),
        _t('CHEBI:61266', 'hemicellulose', 'is a monomer constituent of'),
        _t('CHEBI:171668', '4-O-methyl-alpha-D-glucuronoside',
           'is the free-sugar form of the glycosyl unit in'),
    ]),
    'cellobiose': _m('CHEBI:36217', 'beta-cellobiose', 'exact', [
        _t('CHEBI:17057', 'cellobiose', 'is a'),
        _t('CHEBI:180683', 'lignocellulose',
           'is a structural fragment of the cellulose component of'),
    ]),
    'cellotetraose': _m('CHEBI:62974', 'cellotetraose', 'stereo_unspecified', [
        _t('CHEBI:3523', 'cellodextrin', 'is a'),
    ]),
    'cellohexaose': _m('CHEBI:49533', 'beta-D-cellohexaose', 'exact', [
        _t('CHEBI:3523', 'cellodextrin', 'is a'),
    ]),
    'cellooctaose': _m('CHEBI:49542', 'cellooctaose', 'exact', [
        _t('CHEBI:3523', 'cellodextrin', 'is a'),
    ]),
    'ggm_compact': _m('CHEBI:24144', 'galactoglucomannan', 'broader', [
        _t('CHEBI:53463', 'heptasaccharide', 'is a'),
        _t('CHEBI:61298', 'D-glucose- and/or D-galactose-substituted mannan', 'is a fragment of'),
        _t('GO:0046354', 'mannan biosynthetic process', 'produced by'),
    ]),
    'ggm_extended': _m('CHEBI:24144', 'galactoglucomannan', 'broader', [
        _t('CHEBI:61298', 'D-glucose- and/or D-galactose-substituted mannan', 'is a fragment of'),
        _t('CHEBI:180683', 'lignocellulose', 'is a component of'),
        _t('GO:0046354', 'mannan biosynthetic process', 'produced by'),
    ]),
    'agx_compact': _m('CHEBI:28302', 'glucuronoarabinoxylan', 'broader', [
        _t('CHEBI:60938', 'glucuronoxylan', 'is a fragment of'),
        _t('CHEBI:22603', 'arabinoxylans', 'is a fragment of'),
        _t('CHEBI:37166', 'xylan', 'is a fragment of'),
        _t('GO:0045492', 'xylan biosynthetic process', 'produced by'),
    ]),
    'agx_extended': _m('CHEBI:28302', 'glucuronoarabinoxylan', 'broader', [
        _t('CHEBI:60938', 'glucuronoxylan', 'is a representative fragment of'),
        _t('CHEBI:22603', 'arabinoxylans', 'is a representative fragment of'),
        _t('CHEBI:37166', 'xylan', 'is a fragment of'),
        _t('GO:0045492', 'xylan biosynthetic process', 'produced by'),
    ]),
    'pcoumaryl': _m('CHEBI:64555', 'trans-p-coumaryl alcohol', 'exact', [
        _t('CHEBI:64477', 'monolignol', 'has role'),
        _t('CHEBI:28386', '4-hydroxycinnamyl alcohol', 'E-isomer of'),
    ]),
    'coniferyl': _m('CHEBI:17745', 'coniferol', 'exact', [
        _t('CHEBI:64477', 'monolignol', 'has role'),
        _t('CHEBI:64475', 'guaiacyl lignin', 'monomer of (softwood G unit)'),
        _t('GO:1901063', 'guaiacyl lignin biosynthetic process', 'participant in'),
    ]),
    'sinapyl': _m('CHEBI:64557', 'trans-sinapyl alcohol', 'exact', [
        _t('CHEBI:64477', 'monolignol', 'has role'),
        _t('CHEBI:64476', 'syringal lignin', 'monomer of (angiosperm S unit)'),
    ]),
    'lignin_bo4': _m('CHEBI:53650', 'guaiacylglycerol-beta-guaiacyl ether', 'exact', [
        _t('CHEBI:64475', 'guaiacyl lignin', 'is a representative model of'),
        _t('CHEBI:17745', 'coniferol', 'has monomer unit'),
        _t('GO:1901063', 'guaiacyl lignin biosynthetic process', 'linkage formed in'),
    ]),
    'lignin_b5': _m('CHEBI:91184', 'dehydrodiconiferyl alcohol', 'exact', [
        _t('CHEBI:64475', 'guaiacyl lignin', 'is a representative model of'),
        _t('CHEBI:17745', 'coniferol', 'has monomer unit'),
        _t('GO:1901063', 'guaiacyl lignin biosynthetic process', 'linkage formed in'),
    ]),
    'lignin_bb': _m('CHEBI:40', '(+)-pinoresinol', 'exact', [
        _t('CHEBI:8225', 'pinoresinol', 'is a'),
        _t('CHEBI:25036', 'lignan', 'is a'),
        _t('CHEBI:64475', 'guaiacyl lignin',
           'is a representative model of the beta-beta linkage in'),
        _t('CHEBI:76924', 'plant metabolite', 'has role'),
    ]),
    'lignin_55': _m('CHEBI:64475', 'guaiacyl lignin', 'broader', [
        _t('CHEBI:17745', 'coniferol', 'is an oxidative-coupling dimer of'),
        _t('CHEBI:25497', 'neolignan', 'is a'),
        _t('CHEBI:22888', 'biphenyls', 'is a'),
        _t('GO:1901063', 'guaiacyl lignin biosynthetic process', 'linkage formed in'),
    ]),
    'lignin_trimer': _m('CHEBI:64475', 'guaiacyl lignin', 'broader', [
        _t('CHEBI:91184', 'dehydrodiconiferyl alcohol', 'has substructure'),
        _t('CHEBI:91172', 'guaiacylglycerol beta-coniferyl ether',
           'shares beta-O-4 guaiacylglycerol ether linkage motif'),
        _t('CHEBI:17745', 'coniferol', 'has monomer unit'),
        _t('GO:1901063', 'guaiacyl lignin biosynthetic process', 'linkages formed in'),
    ]),
    'o2': _m('CHEBI:15379', 'dioxygen', 'exact', [
        _t('CHEBI:27140', 'triplet dioxygen',
           'is a closed-shell depiction of the ground-state spin form'),
        _t('CHEBI:63248', 'oxidising agent', 'has role'),
    ]),
    'co2': _m('CHEBI:16526', 'carbon dioxide', 'exact', [
        _t('ENVO:01000451', 'atmospheric carbon dioxide', 'is a molecular model of'),
    ]),
    'h2o': _m('CHEBI:15377', 'water', 'exact', [
        _t('ENVO:01000266', 'water vapour', 'is a molecular model of'),
        _t('ENVO:00002006', 'liquid water', 'is a molecular model of'),
    ]),
    'oleic_acid': _m('CHEBI:16196', 'oleic acid', 'exact', [
        _t('CHEBI:25413', 'monounsaturated fatty acid', 'is a'),
        _t('CHEBI:15904', 'long-chain fatty acid', 'is a'),
        _t('CHEBI:76924', 'plant metabolite', 'has role'),
    ]),
    'palmitic_acid': _m('CHEBI:15756', 'hexadecanoic acid', 'exact', [
        _t('CHEBI:39418', 'straight-chain saturated fatty acid', 'is a'),
        _t('CHEBI:26607', 'saturated fatty acid', 'is a representative model of'),
        _t('CHEBI:15904', 'long-chain fatty acid', 'is a'),
    ]),
    'triolein': _m('CHEBI:53753', 'triolein', 'exact', [
        _t('CHEBI:17855', 'triglyceride', 'is a representative model of'),
        _t('GO:0019432', 'triglyceride biosynthetic process', 'produced by'),
    ]),
    'tripalmitin': _m('CHEBI:77393', 'tripalmitin', 'exact', [
        _t('CHEBI:17855', 'triglyceride', 'is a representative model of'),
        _t('CHEBI:15756', 'hexadecanoic acid', 'has acyl component'),
        _t('GO:0019432', 'triglyceride biosynthetic process', 'produced by'),
    ]),
    'alpha_pinene': _m('CHEBI:28660', '(-)-alpha-pinene', 'exact', [
        _t('CHEBI:36740', 'alpha-pinene', 'is a'),
        _t('CHEBI:35187', 'monoterpene', 'is a'),
        _t('CHEBI:17187', 'pinene', 'is a'),
        _t('CHEBI:76924', 'plant metabolite', 'has role'),
        _t('GO:0016099', 'monoterpenoid biosynthetic process', 'produced by'),
    ]),
    'beta_pinene': _m('CHEBI:28359', '(-)-beta-pinene', 'exact', [
        _t('CHEBI:50025', 'beta-pinene', 'is a'),
        _t('CHEBI:35187', 'monoterpene', 'is a'),
        _t('ENVO:03510018', 'resin', 'is a volatile component of'),
        _t('GO:0016099', 'monoterpenoid biosynthetic process', 'produced by'),
    ]),
    'abietic_acid': _m('CHEBI:28987', 'abietic acid', 'exact', [
        _t('CHEBI:36762', 'abietane diterpenoid', 'is a'),
        _t('CHEBI:753451', 'rosin', 'is the principal constituent of'),
        _t('CHEBI:76924', 'plant metabolite', 'has role'),
        _t('GO:0016102', 'diterpenoid biosynthetic process', 'produced by'),
    ]),
    'charcoal': _m('CHEBI:29863', 'coronene', 'exact', [
        _t('ENVO:01000560', 'charcoal', 'is a representative model of'),
        _t('CHEBI:91090', 'charcoal', 'is a representative model of'),
        _t('CHEBI:35300', 'ortho- and peri-fused polycyclic arene', 'is a'),
    ]),
    'lignite': _m('CHEBI:25238', 'methoxybenzoic acid', 'broader', [
        _t('ENVO:2000008', 'lignite', 'is a representative model of'),
        _t('ENVO:02000091', 'coal', 'is a representative model fragment of'),
        _t('CHEBI:64475', 'guaiacyl lignin', 'is derived from'),
        _t('CHEBI:67649', '(+)-(7S,8S)-guaiacylglycerol-beta-vanillic acid ether',
           'is a close structural analogue of'),
    ]),
    'sub_bituminous': _m('CHEBI:59659', 'pyrenes', 'broader', [
        _t('ENVO:2000009', 'sub-bitumous coal', 'is a representative model of'),
        _t('ENVO:02000091', 'coal', 'is a model fragment of'),
        _t('CHEBI:34093', '1-Hydroxypyrene', 'is a 6-ethyl derivative of'),
    ]),
    'bituminous': _m('CHEBI:33853', 'phenols', 'broader', [
        _t('ENVO:2000010', 'bitumous coal', 'is a representative model of'),
        _t('ENVO:02000091', 'coal', 'is a representative model of'),
        _t('CHEBI:51348', 'tetraphene', 'has part'),
        _t('CHEBI:35230', 'fossil fuel', 'has role'),
    ]),
    'anthracite': _m('CHEBI:35300', 'ortho- and peri-fused polycyclic arene', 'broader', [
        _t('ENVO:2000011', 'anthracite', 'is a representative model of'),
        _t('ENVO:02000091', 'coal', 'is a fragment of'),
        _t('CHEBI:33418', 'graphite', 'is a defect-containing finite-sheet approximation of'),
    ]),
    'kerogen': _m('CHEBI:38836', '1-benzothiophenes', 'broader', [
        _t('ENVO:01001231', 'kerogen', 'is a representative model fragment of'),
        _t('ENVO:01001235', 'planktonic kerogen', 'is a representative model of'),
        _t('CHEBI:38180', 'polycyclic heteroarene', 'is a'),
        _t('CHEBI:76224', 'aromatic ketone', 'is a'),
    ]),
    'crude_oil': _m('CHEBI:35300', 'ortho- and peri-fused polycyclic arene', 'broader', [
        _t('ENVO:00002984', 'petroleum', 'representative model of (asphaltene)'),
        _t('CHEBI:29863', 'coronene', 'aromatic core of the model'),
        _t('CHEBI:35230', 'fossil fuel', 'has role'),
    ]),
    'udp_glucose': _m('CHEBI:46229', 'UDP-alpha-D-glucose', 'exact', [
        _t('CHEBI:18066', 'UDP-D-glucose', 'is a representative model of'),
        _t('CHEBI:61109', 'pyrimidine nucleotide-sugar', 'is a'),
        _t('CHEBI:18246', '(1->4)-beta-D-glucan', 'is the activated glucosyl donor of'),
        _t('GO:0030244', 'cellulose biosynthetic process', 'participates in'),
        _t('GO:0009225', 'nucleotide-sugar metabolic process', 'participant in'),
    ]),
    'gdp_mannose': _m('CHEBI:15820', 'GDP-alpha-D-mannose', 'exact', [
        _t('CHEBI:21168', 'GDP-mannose', 'is a stereo-unspecified model of'),
        _t('CHEBI:64711', 'purine nucleotide-sugar', 'is a'),
        _t('CHEBI:24144', 'galactoglucomannan',
           'is the activated mannosyl donor for biosynthesis of'),
        _t('GO:0009225', 'nucleotide-sugar metabolic process', 'participant in'),
        _t('GO:0046354', 'mannan biosynthetic process', 'glycosyl donor for'),
    ]),
    'phenylalanine': _m('CHEBI:17295', 'L-phenylalanine', 'exact', [
        _t('CHEBI:33856', 'aromatic amino acid', 'is a'),
        _t('CHEBI:83813', 'proteinogenic amino acid', 'is a'),
        _t('CHEBI:26004', 'phenylpropanoid', 'is the biosynthetic precursor of'),
        _t('GO:0009699', 'phenylpropanoid biosynthetic process', 'precursor in'),
    ]),
    'acetyl_coa': _m('CHEBI:15351', 'acetyl-CoA', 'exact', [
        _t('CHEBI:17984', 'acyl-CoA', 'is a'),
        _t('GO:0006633', 'fatty acid biosynthetic process', 'is a substrate/building block of'),
        _t('GO:0019287', 'isopentenyl diphosphate biosynthetic process, mevalonate pathway',
           'is a precursor in'),
    ]),
    'ipp': _m('CHEBI:16584', 'isopentenyl diphosphate', 'exact', [
        _t('CHEBI:37841', 'isoprenoid phosphate', 'is a member of the class'),
        _t('CHEBI:24913', 'isoprenoid', 'is the C5 building-block precursor of'),
        _t('GO:0009240', 'isopentenyl diphosphate biosynthetic process', 'is the output of'),
        _t('CHEBI:76924', 'plant metabolite', 'has role'),
    ]),
    'dmapp': _m('CHEBI:16057', 'prenyl diphosphate', 'exact', [
        _t('CHEBI:26250', 'prenol phosphate', 'is a'),
        _t('CHEBI:37841', 'isoprenoid phosphate', 'is a'),
        _t('GO:0050992', 'dimethylallyl diphosphate biosynthetic process', 'product of'),
    ]),
}


# --------------------------------------------------------------------------------------
def terms_for(mol_id: str) -> dict:
    """All annotation layers for one molecule, each term carrying its IRI."""
    rec = molecules.BY_ID[mol_id]
    m = MOLECULE_TERMS[mol_id]
    add = lambda t: {**t, "iri": iri(t["curie"])}
    return {
        "structure": {**add(m["structure"]), "match": m["structure"]["match"]},
        "represents": [add(t) for t in m.get("represents", [])],
        "category": [add(t) for t in CATEGORY_TERMS.get(rec["category"], [])],
        "project": [add(t) for t in PROJECT_TERMS],
    }


def _wrap(text: str, width: int) -> list[str]:
    """Greedy word wrap; a single word longer than `width` is hard-split."""
    words = [w[i:i + width] for w in text.split() for i in range(0, len(w), width)]
    lines, cur = [], ""
    for w in words:
        if cur and len(cur) + 1 + len(w) > width:
            lines.append(cur)
            cur = w
        else:
            cur = f"{cur} {w}" if cur else w
    return lines + [cur] if cur else lines


# Fixed-column layout (cols 1-11 "REMARK 950 ", 12-22 tag, 23-39 CURIE, 40-80 text):
#   REMARK 950 STRUCTURE  CHEBI:49533      beta-D-cellohexaose [match=exact]
#   REMARK 950 CATEGORY   GO:0052324       plant-type cell wall cellulose
#   REMARK 950                             biosynthetic process [produced by]
# A line with a tag starts a term; a blank tag/CURIE column continues the previous one.
# IRI = http://purl.obolibrary.org/obo/<PREFIX>_<ID> of the CURIE (see iri()).
_TAG_W, _CURIE_W = 11, 17
_TEXT_W = PAYLOAD - _TAG_W - _CURIE_W


def _term_lines(tag: str, t: dict) -> list[str]:
    rel = f" [{t['relation']}]" if t.get("relation") else ""
    text = _wrap(f"{t['label']}{rel}", _TEXT_W)
    first = f"{tag:<{_TAG_W}}{t['curie']:<{_CURIE_W}}{text[0]}"
    return [first] + [" " * (_TAG_W + _CURIE_W) + ln for ln in text[1:]]


def remark_block(mol_id: str) -> list[str]:
    """The COMPND + REMARK 950 header lines (each <= 80 cols, no trailing newline)."""
    rec = molecules.BY_ID[mol_id]
    a = terms_for(mol_id)
    s = a["structure"]
    body = ["OBO FOUNDRY ONTOLOGY ANNOTATIONS (energetics-pymol ontology.py)",
            "resolved via EMBL-EBI OLS4; IRI = http://purl.obolibrary.org/obo/",
            "  + CURIE with ':' as '_'. STRUCTURE match levels: exact,",
            "  stereo_unspecified, close, broader (see ontology.MATCH_LEVELS)"]
    body += _term_lines("STRUCTURE", {**s, "relation": f"match={s['match']}"})
    for t in a["represents"]:
        body += _term_lines("REPRESENTS", t)
    for t in a["category"]:
        body += _term_lines("CATEGORY", t)
    for t in a["project"]:
        body += _term_lines("PROJECT", t)
    # wwPDB COMPND tokens (MOL_ID / MOLECULE). RDKit reads the name into _Name; Mol* ignores
    # COMPND without a CHAIN token (these are single HETATM residues), so it is informational.
    # PDB is ASCII: transliterate dashes; ';' and ':' are COMPND token delimiters.
    name = rec["name"].replace("\u2014", "-").replace("\u2013", "-").replace(";", ",")
    name = name.replace(":", "-")
    name = name.encode("ascii", "ignore").decode()
    mol = _wrap(f"MOLECULE: {name};", 69)
    compnd = ["COMPND    MOL_ID: 1;"] + [f"COMPND {i + 2:>3d} {ln}" for i, ln in enumerate(mol)]
    lines = compnd + [REMARK + ln for ln in body]
    assert all(ln.isascii() for ln in lines), f"{mol_id}: non-ASCII header"
    too_long = [ln for ln in lines if len(ln) > 80]
    assert not too_long, f"{mol_id}: PDB lines > 80 cols: {too_long}"
    return lines


def annotate_pdb(path: str, mol_id: str) -> bool:
    """Idempotently (re)write the COMPND/REMARK 950 header. Returns True if the file changed."""
    with open(path) as fh:
        old = fh.read().splitlines()
    kept = [ln for ln in old if not (ln.startswith("COMPND") or ln.startswith("REMARK 950"))]
    new = remark_block(mol_id) + kept
    if new == old:
        return False
    with open(path, "w") as fh:
        fh.write("\n".join(new) + "\n")
    return True


def check() -> list[str]:
    errs = []
    for rec in molecules.MANIFEST:
        mid = rec["id"]
        if mid not in MOLECULE_TERMS:
            errs.append(f"{mid}: no MOLECULE_TERMS entry")
            continue
        if MOLECULE_TERMS[mid]["structure"]["match"] not in MATCH_LEVELS:
            errs.append(f"{mid}: bad match level")
        if rec["category"] not in CATEGORY_TERMS:
            errs.append(f"{mid}: category {rec['category']!r} has no CATEGORY_TERMS")
        path = os.path.join(PDB_DIR, f"{mid}.pdb")
        if os.path.exists(path):
            with open(path) as fh:
                lines = fh.read().splitlines()
            if [ln for ln in lines if ln.startswith(("COMPND", "REMARK 950"))] != remark_block(mid):
                errs.append(f"{mid}: pdb REMARK 950 block missing or stale (run --annotate)")
            if any(len(ln) > 80 for ln in lines):
                errs.append(f"{mid}: pdb line > 80 columns")
    extra = set(MOLECULE_TERMS) - set(molecules.BY_ID)
    errs += [f"{x}: annotated but not in manifest" for x in sorted(extra)]
    return errs


def _ols_term(curie: str) -> dict:
    """Fetch one term from the EMBL-EBI OLS4 REST API (stdlib only)."""
    import json
    import urllib.parse
    import urllib.request
    ont = curie.split(":", 1)[0].lower()
    enc = urllib.parse.quote(urllib.parse.quote(iri(curie), safe=""), safe="")
    url = f"https://www.ebi.ac.uk/ols4/api/ontologies/{ont}/terms/{enc}"
    with urllib.request.urlopen(url, timeout=30) as fh:
        r = json.load(fh)
    ik = (r.get("annotation", {}).get("inchi_key_string") or [None])[0]
    return dict(label=r.get("label"), obsolete=bool(r.get("is_obsolete")), inchikey=ik)


def verify_ols() -> list[str]:
    """Network gate: every CURIE exists in OLS4 with exactly our label and is not obsolete;
    every structure term's InChIKey agrees with the model at its declared match level."""
    from concurrent.futures import ThreadPoolExecutor
    from rdkit import Chem
    import build_pdb

    terms = {t["curie"]: t["label"] for t in PROJECT_TERMS}
    terms.update({t["curie"]: t["label"] for ts in CATEGORY_TERMS.values() for t in ts})
    for m in MOLECULE_TERMS.values():
        terms[m["structure"]["curie"]] = m["structure"]["label"]
        terms.update({t["curie"]: t["label"] for t in m["represents"]})
    with ThreadPoolExecutor(12) as ex:
        live = dict(zip(terms, ex.map(_ols_term, terms)))
    errs = [f"{c}: OLS label {live[c]['label']!r} != {lab!r}" for c, lab in terms.items()
            if live[c]["label"] != lab]
    errs += [f"{c}: obsolete in OLS" for c in terms if live[c]["obsolete"]]
    for rec in molecules.MANIFEST:
        s = MOLECULE_TERMS[rec["id"]]["structure"]
        if s["match"] not in ("exact", "stereo_unspecified", "close"):
            continue
        model = Chem.MolToInchiKey(build_pdb.get_mol(rec))
        term = live[s["curie"]]["inchikey"]
        if not term:
            errs.append(f"{rec['id']}: {s['curie']} has no InChIKey to check a {s['match']} match")
        elif s["match"] == "exact" and term != model:
            errs.append(f"{rec['id']}: exact match but InChIKey {term} != model {model}")
        elif s["match"] != "exact" and term[:14] != model[:14]:
            errs.append(f"{rec['id']}: {s['match']} match but constitution differs ({term} vs {model})")
    print(f"OLS4: {len(terms)} terms checked")
    return errs


def main(argv):
    if "--verify-ols" in argv:
        errs = verify_ols()
        for e in errs:
            print("[FAIL]", e)
        print(f"OLS verification: {len(errs)} problem(s)")
        return 1 if errs else 0
    if "--annotate" in argv:
        n = sum(annotate_pdb(os.path.join(PDB_DIR, f"{r['id']}.pdb"), r["id"])
                for r in molecules.MANIFEST if os.path.exists(os.path.join(PDB_DIR, f"{r['id']}.pdb")))
        print(f"annotated {n} pdb file(s)")
    if "--check" in argv or "--annotate" in argv:
        errs = check()
        for e in errs:
            print("[FAIL]", e)
        print(f"ontology check: {len(errs)} problem(s)")
        return 1 if errs else 0
    for rec in molecules.MANIFEST:
        s = MOLECULE_TERMS.get(rec["id"], {}).get("structure")
        print(f"{rec['id']:16s} {s['curie'] if s else '-':14s} {s['match'] if s else '':20s} "
              f"{s['label'] if s else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
