"""
molecules.py — verified molecule manifest for the "Ecosystems as energy fields"
(Falk, Swetnam & McKenzie 2026) figure set.

This is the SINGLE SOURCE OF TRUTH. Every structure here was cross-checked against
PubChem / ChEBI by an adversarial verification pass. Where chemistry differs by clade
we use SOFTWOOD / CONIFER (Pinus) representatives, because the study system is
Pinus ponderosa.

Each manifest record is a dict:
    id        unique slug -> output file stem (pdb/<id>.pdb, renders/<id>.png)
    name      human-readable name (figure caption)
    resname   3-char PDB residue code (CCD style where one exists)
    category  grouping for the catalog/README
    formula   expected molecular formula (Hill); build_pdb asserts this
    smiles    isomeric SMILES, OR None if `builder` supplies an RDKit Mol
    builder   name of a function in polymers.py (string) that returns a Mol, OR None
    size      'na' | 'repeat' | 'compact' | 'extended'
    show_h    keep explicit hydrogens visible in the default render
    note      provenance / caption note

`smiles` strings are always run through RDKit with isomeric stereo preserved.
Polymers (cellulose oligomers, softwood hemicellulose, lignin oligomers) are built in
polymers.py via fragment-then-join so the inter-unit CONECT bonds are real bonds.
"""

# --------------------------------------------------------------------------------------
# Verified monomer SMILES reused by polymers.py (kept here so there is ONE definition).
# polymers.py asserts its atom-mapped building blocks canonicalize to exactly these.
# --------------------------------------------------------------------------------------
MONOMER_SMILES = {
    "bdglc":  "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",        # beta-D-glucopyranose  C6H12O6
    "bdman":  "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",       # beta-D-mannopyranose  C6H12O6
    "bdgal":  "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",         # beta-D-galactopyranose C6H12O6
    "adgal":  "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",          # alpha-D-galactopyranose (branch donor)
    "bdxyl":  "C1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O",             # beta-D-xylopyranose   C5H10O5
    "alaraf": "C([C@H]1[C@@H]([C@H]([C@@H](O1)O)O)O)O",             # alpha-L-arabinofuranose C5H10O5
    "meglca": "CO[C@H]1[C@@H]([C@H]([C@H](O[C@@H]1C(=O)O)O)O)O",    # 4-O-methyl-alpha-D-glucuronic acid C7H12O7
}


def _m(id, name, resname, category, formula, smiles=None, builder=None,
       size="na", show_h=False, note=""):
    return dict(id=id, name=name, resname=resname, category=category, formula=formula,
                smiles=smiles, builder=builder, size=size, show_h=show_h, note=note)


MANIFEST = [

    # ---- Sugar / polysaccharide monomers -------------------------------------------
    _m("bdglc", "beta-D-glucopyranose", "BGC", "Sugar monomers", "C6H12O6",
       smiles=MONOMER_SMILES["bdglc"],
       note="Cellulose building block & primary photosynthate (PubChem CID 64689)."),
    _m("glc_open", "D-glucose (open-chain, Fischer aldehyde)", "GLC", "Sugar monomers", "C6H12O6",
       smiles="OC[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C=O",
       note="Acyclic aldohexose: the photosynthesis product before ring closure (ChEBI aldehydo-D-glucose)."),
    _m("bdman", "beta-D-mannopyranose", "MAN", "Sugar monomers", "C6H12O6",
       smiles=MONOMER_SMILES["bdman"],
       note="Principal backbone sugar of softwood galactoglucomannan (C2 epimer of glucose)."),
    _m("bdgal", "beta-D-galactopyranose", "GAL", "Sugar monomers", "C6H12O6",
       smiles=MONOMER_SMILES["bdgal"],
       note="Galactoglucomannan side-branch sugar (C4 epimer of glucose)."),
    _m("bdxyl", "beta-D-xylopyranose", "XYS", "Sugar monomers", "C5H10O5",
       smiles=MONOMER_SMILES["bdxyl"],
       note="Pentose backbone of softwood arabinoglucuronoxylan (no C6 CH2OH)."),
    _m("alaraf", "alpha-L-arabinofuranose", "ARA", "Sugar monomers", "C5H10O5",
       smiles=MONOMER_SMILES["alaraf"],
       note="Furanose (5-membered) side branch of arabinoglucuronoxylan (PubChem CID 7044039)."),
    _m("meglca", "4-O-methyl-alpha-D-glucuronic acid", "MGA", "Sugar monomers", "C7H12O7",
       smiles=MONOMER_SMILES["meglca"],
       note="Acidic side branch of softwood xylan; C6 oxidized to COOH, O4 methylated (PubChem CID 445994)."),

    # ---- Cellulose: monomer -> repeat -> oligomers (compact + extended) -------------
    _m("cellobiose", "Cellobiose (beta-D-Glcp-(1->4)-D-Glcp)", "BGC", "Cellulose", "C12H22O11",
       builder="cellulose:2", size="repeat",
       note="The cellulose repeating disaccharide; shows the beta-1,4 glycosidic bond."),
    _m("cellotetraose", "Cellotetraose (cellulose, N=4)", "BGC", "Cellulose", "C24H42O21",
       builder="cellulose:4", size="compact",
       note="Compact cellulose oligomer."),
    _m("cellohexaose", "Cellohexaose (cellulose, N=6)", "BGC", "Cellulose", "C36H62O31",
       builder="cellulose:6", size="compact",
       note="Standard crystalline-cellulose surrogate; twofold-screw ribbon w/ O3-H...O5 H-bonds."),
    _m("cellooctaose", "Cellooctaose (cellulose, N=8)", "BGC", "Cellulose", "C48H82O41",
       builder="cellulose:8", size="extended",
       note="Extended cellulose chain for a longer-ribbon figure."),

    # ---- Hemicellulose: softwood galactoglucomannan & arabinoglucuronoxylan ---------
    _m("ggm_compact", "Galactoglucomannan fragment (compact)", "GGM", "Hemicellulose", None,
       builder="ggm:compact", size="compact",
       note="beta-1,4 Man/Glc backbone + alpha-1,6 Gal branch; principal softwood hemicellulose."),
    _m("ggm_extended", "Galactoglucomannan fragment (extended, acetylated)", "GGM", "Hemicellulose", None,
       builder="ggm:extended", size="extended",
       note="Longer GGM backbone with two Gal branches and native O-acetyl groups."),
    _m("agx_compact", "Arabinoglucuronoxylan fragment (compact)", "AGX", "Hemicellulose", None,
       builder="agx:compact", size="compact",
       note="beta-1,4 xylan backbone + alpha-1,2 MeGlcA + alpha-1,3 arabinofuranose branches."),
    _m("agx_extended", "Arabinoglucuronoxylan fragment (extended)", "AGX", "Hemicellulose", None,
       builder="agx:extended", size="extended",
       note="Longer xylan backbone with multiple acidic and arabinose branches."),

    # ---- Lignin: monolignols + softwood (guaiacyl) linkage models -------------------
    _m("pcoumaryl", "p-Coumaryl alcohol (H monolignol)", "LIG", "Lignin", "C9H10O2",
       smiles="OC1=CC=C(/C=C/CO)C=C1",
       note="Trace H units in softwood; included for H/G/S contrast (PubChem CID 5280535)."),
    _m("coniferyl", "Coniferyl alcohol (G monolignol)", "LIG", "Lignin", "C10H12O3",
       smiles="COC1=CC(=CC=C1O)/C=C/CO",
       note="THE primary softwood (guaiacyl) monolignol (PubChem CID 1549095)."),
    _m("sinapyl", "Sinapyl alcohol (S monolignol)", "LIG", "Lignin", "C11H14O4",
       smiles="COC1=CC(=CC(=C1O)OC)/C=C/CO",
       note="Syringyl monomer; essentially absent from softwood -- caption as the angiosperm contrast (CID 5280507)."),
    _m("lignin_bo4", "Guaiacylglycerol-beta-guaiacyl ether (beta-O-4)", "LIG", "Lignin", "C17H20O6",
       smiles="COC1=CC=CC=C1OC(CO)C(O)C2=CC=C(O)C(OC)=C2",
       note="The dominant softwood interunit linkage (~50%); classic beta-O-4 model (PubChem CID 6424189)."),
    _m("lignin_b5", "Dehydrodiconiferyl alcohol (beta-5 phenylcoumaran)", "LIG", "Lignin", "C20H22O6",
       smiles="COc1cc(/C=C/CO)cc2c1OC(c1ccc(O)c(OC)c1)C2CO",
       note="beta-5 phenylcoumaran dimer, both rings guaiacyl (PubChem CID 5372367)."),
    _m("lignin_bb", "Pinoresinol (beta-beta resinol)", "LIG", "Lignin", "C20H22O6",
       smiles="COC1=C(C=CC(=C1)[C@@H]2[C@H]3CO[C@@H]([C@H]3CO2)C4=CC(=C(C=C4)O)OC)O",
       show_h=False,
       note="beta-beta furofuran dimer; also a genuine conifer lignan (keep defined stereo)."),
    _m("lignin_55", "5-5' biphenyl diconiferyl alcohol", "LIG", "Lignin", "C20H22O6",
       smiles="OC/C=C/c1cc(OC)c(O)c(-c2cc(/C=C/CO)cc(OC)c2O)c1",
       size="extended",
       note="5-5 biphenyl C-C linkage (most common softwood C-C bond, often within dibenzodioxocin)."),
    _m("lignin_trimer", "G-(beta-O-4)-G-(beta-5)-G lignin trimer", "LIG", "Lignin", None,
       builder="lignin:trimer", size="extended",
       note="Defined representative guaiacyl trimer combining beta-O-4 + beta-5 (caption as a model, not bulk lignin)."),

    # ---- Combustion small molecules (the HHV / enthalpy-of-combustion framing) ------
    _m("o2", "Dioxygen (O2)", "OXY", "Combustion", "O2",
       smiles="O=O", show_h=False,
       note="Oxidant in biomass + O2 -> CO2 + H2O (drawn closed-shell; ground state is a triplet diradical)."),
    _m("co2", "Carbon dioxide (CO2)", "CO2", "Combustion", "CO2",
       smiles="O=C=O", show_h=False,
       note="Combustion product / photosynthesis substrate."),
    _m("h2o", "Water (H2O)", "HOH", "Combustion", "H2O",
       smiles="O", show_h=True,
       note="Combustion product; the HHV/LHV distinction hinges on its heat of vaporization."),

    # ---- Fats (energy-dense biomass constituent named in the paper) -----------------
    _m("oleic_acid", "Oleic acid (18:1 cis-9)", "OLA", "Fats", "C18H34O2",
       smiles=r"CCCCCCCC/C=C\CCCCCCCC(=O)O",
       note="Major fatty acid of pine oleoresin/tall oil; cis double bond gives the chain kink."),
    _m("palmitic_acid", "Palmitic acid (16:0)", "PLM", "Fats", "C16H32O2",
       smiles="CCCCCCCCCCCCCCCC(=O)O",
       note="Common saturated fatty acid; straight-chain contrast to oleic."),
    _m("triolein", "Triolein (glyceryl trioleate)", "TRO", "Fats", "C57H104O6",
       smiles=r"CCCCCCCC/C=C\CCCCCCCC(=O)OCC(OC(=O)CCCCCCC/C=C\CCCCCCCC)COC(=O)CCCCCCC/C=C\CCCCCCCC",
       note="Representative triacylglyceride (the storage form of 'fats'); three oleate tails on glycerol."),
    _m("tripalmitin", "Tripalmitin (glyceryl tripalmitate)", "TPM", "Fats", "C51H98O6",
       smiles="CCCCCCCCCCCCCCCC(=O)OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC",
       note="Saturated triacylglyceride alternative."),

    # ---- Extractives (conifer; paper notes extractives at 10-45%) -------------------
    _m("alpha_pinene", "alpha-Pinene", "APN", "Extractives", "C10H16",
       smiles="CC1=CC[C@H]2C[C@@H]1C2(C)C",
       note="Dominant Pinus monoterpene; (1S,5S)-(-) enantiomer common in N. American pines (CID 440968)."),
    _m("beta_pinene", "beta-Pinene", "BPN", "Extractives", "C10H16",
       smiles="CC1([C@H]2CCC(=C)[C@@H]1C2)C",
       note="Exocyclic-methylene isomer of alpha-pinene; co-occurs in oleoresin (CID 440967)."),
    _m("abietic_acid", "Abietic acid", "ABT", "Extractives", "C20H30O2",
       smiles="CC(C)C1=CC2=CC[C@@H]3[C@@]([C@H]2CC1)(CCC[C@@]3(C)C(=O)O)C",
       note="Signature conifer diterpene resin acid (rosin); conjugated 7,13-diene (PubChem CID 10569)."),

    # ---- Fossil carbon: representative MODELS for the long-term carbon vaults ----------
    # Coal, char, kerogen and oil are heterogeneous macromolecular solids/mixtures, not
    # single compounds. Each entry is a defensible representative MODEL (verified, RDKit-
    # buildable). Across coal rank, aromaticity rises and O falls: lignite -> anthracite.
    _m("charcoal", "Charcoal — polycyclic aromatic cluster", "CHA", "Fossil carbon", "C24H12",
       smiles="c1cc2ccc3ccc4ccc5ccc6ccc1c1c2c3c4c5c61",
       note="Representative model: coronene, a small graphene-like fused-aromatic sheet. Real charcoal/"
            "pyrogenic carbon is defect-rich, polydisperse fused aromatic clusters with residual O and ash."),
    _m("lignite", "Lignite — O-rich coal fragment", "LGN", "Fossil carbon", "C17H18O7",
       smiles="COc1cc(C(=O)O)ccc1OCC(O)c1ccc(O)c(OC)c1",
       note="Representative model: lowest-rank coal, still lignin-derived — small aromatics with methoxyl, "
            "phenol, carboxyl and ether/aliphatic bridges (high O/C). The O-rich end of the rank series."),
    _m("sub_bituminous", "Sub-bituminous — intermediate coal fragment", "SBC", "Fossil carbon", "C18H14O",
       smiles="CCc1ccc2ccc3c(O)ccc4ccc1c2c34",
       note="Representative model: a hydroxy-pyrene with a short alkyl chain — larger aromatic unit and far "
            "less oxygen than lignite, sitting between lignite and bituminous in rank."),
    _m("bituminous", "Bituminous — classic coal model fragment", "BIT", "Fossil carbon", "C29H24O",
       smiles="CC1Cc2ccc3ccc4cc5ccccc5cc4c3c2C1Cc1ccc(O)cc1",
       note="Representative model (Shinn/Wiser-style): fused aromatic clusters linked by short methylene "
            "bridges with a hydroaromatic ring and trace oxygen — the hallmark mid-rank coal motif."),
    _m("anthracite", "Anthracite — near-graphitic aromatic sheet", "ANT", "Fossil carbon", "C38H16",
       smiles="c1cc2ccc3ccc4ccc5ccc6ccc7ccc8ccc1c1c2c2c3c4c3c5c6c4c7c8c1c2c34",
       note="Representative model: a large peri-fused PAH (circumcoronene class) — highest-rank coal, "
            "near-pure fused aromatic carbon approaching graphite (minimal H, ~no O)."),
    _m("kerogen", "Kerogen — Type II (marine) fragment", "KER", "Fossil carbon", "C26H28OS",
       smiles="CCCCCc1cc2ccc3ccc4sc(C(=O)CCCC)cc4c3c2cc1",
       note="Representative model: a fused-aromatic core with aliphatic chains, a thiophene (organic sulfur) "
            "and a carbonyl — the insoluble organic matter that cracks to petroleum on burial."),
    _m("crude_oil", "Crude oil — island asphaltene archetype", "OIL", "Fossil carbon", "C42H48",
       smiles="CCCCCCc1cc2ccc3c(CCCCCC)cc4ccc5c(CCCCCC)cc6ccc1c1c6c5c4c3c21",
       note="Representative model: a 'continental/island' asphaltene (Yen-Mullins) — a polycyclic aromatic "
            "core with alkyl chains, the canonical heavy molecule of crude oil (a complex mixture)."),
]

# id -> record, for convenience
BY_ID = {m["id"]: m for m in MANIFEST}

CATEGORIES = ["Sugar monomers", "Cellulose", "Hemicellulose", "Lignin",
              "Combustion", "Fats", "Extractives", "Fossil carbon"]

if __name__ == "__main__":
    # quick textual catalog (no RDKit needed)
    for cat in CATEGORIES:
        print(f"\n## {cat}")
        for m in MANIFEST:
            if m["category"] == cat:
                src = m["smiles"] if m["smiles"] else f"<builder {m['builder']}>"
                print(f"  {m['id']:16s} {m['formula'] or '?':10s} {m['name']}")
    print(f"\nTotal: {len(MANIFEST)} molecules")
