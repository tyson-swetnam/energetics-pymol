# alpha-Pinene
# Dominant Pinus monoterpene; (1S,5S)-(-) enantiomer common in N. American pines (CID 440968).
reinitialize
bg_color white
load /Users/tswetnam/Desktop/energetics-pymol/pdb/alpha_pinene.pdb, mol
set connect_mode, 1          # trust CONECT records (non-standard residues)

hide everything
show sticks, mol
hide everything, (mol and hydro)
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
zoom mol, 2.0
ray 1600, 1200
png /Users/tswetnam/Desktop/energetics-pymol/renders/alpha_pinene.png, dpi=300, ray=1
