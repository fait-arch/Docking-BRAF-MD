# Figuras 3D del reporte (PyMOL). Ejecutar desde scripts/:  pymol -cq figuras_pymol.pml
load ../md_prueba/04_prod_004.pdb, cx
remove not (polymer.protein or resn LIG)
alter polymer.protein, resi=str(int(resi)+447)
sort
set ray_opaque_background, 1
bg_color white
set ray_shadows, 0
set antialias, 2
set ray_trace_mode, 1
set ray_trace_color, 0x444444
set cartoon_fancy_helices, 1
set cartoon_side_chain_helper, 1
set valence, 0
set_color c_gris, [0.84, 0.84, 0.82]
set_color c_azul, [0.165, 0.471, 0.839]
set_color c_naranja, [0.922, 0.408, 0.204]
set_color c_aqua, [0.106, 0.686, 0.478]
set_color c_amar, [0.929, 0.631, 0.0]
set_color c_mag, [0.910, 0.482, 0.643]
set_color c_viol, [0.384, 0.314, 0.839]
hide everything
show cartoon, polymer.protein
color c_gris, polymer.protein
color c_amar, polymer.protein and resi 464-471
color c_viol, polymer.protein and resi 492-507
color c_aqua, polymer.protein and resi 529-535
color c_mag, polymer.protein and resi 594-623
show spheres, resn LIG and not hydro
color c_naranja, resn LIG and elem C
color red, resn LIG and elem O
color blue, resn LIG and elem N
orient polymer.protein
turn z, -90
zoom polymer.protein, -3
ray 1500, 1500
png ../resultados_prueba/fig_complejo_general.png, dpi=300

# ------------------------------------------------ acercamiento
hide spheres
set cartoon_color, c_gris
set cartoon_transparency, 0.72
set stick_radius, 0.17
select bolsillo, polymer.protein and resi 471+483+514+529+531+532+583+594 and not name N+C+O
show sticks, bolsillo and not hydro
show sticks, polymer.protein and resi 532 and name N+CA+C+O
color c_azul, (bolsillo or (polymer.protein and resi 532)) and elem C
color red, bolsillo and elem O
color blue, bolsillo and elem N
color yellow, bolsillo and elem S
show sticks, resn LIG and not hydro
show sticks, resn LIG and hydro and neighbor elem O
color white, resn LIG and hydro
set stick_radius, 0.27, resn LIG
distance hb, resn LIG and name O1, polymer.protein and resi 529 and name OG1, 3.5
hide labels, hb
set dash_color, 0x111111
set dash_radius, 0.07
set dash_gap, 0.22
set label_size, -1.55
set float_labels, 1
set label_outline_color, white
set ray_label_specular, 0
set label_color, 0x111111
set label_outline_color, white
set label_font_id, 7
set label_position, (0, 0, 6)
label polymer.protein and resi 483 and name CE, "Lys483"
label polymer.protein and resi 529 and name CG2, "Thr529"
label polymer.protein and resi 532 and name SG, "Cys532"
label polymer.protein and resi 531 and name CE3, "Trp531"
label polymer.protein and resi 594 and name OD2, "Asp594"
label polymer.protein and resi 583 and name CZ, "Phe583"
label polymer.protein and resi 514 and name CD2, "Leu514"
label polymer.protein and resi 471 and name CG1, "Val471"
set label_position, (-2.6, -0.6, 6), polymer.protein and resi 529 and name CG2
set label_position, (-1.8, -1.6, 6), polymer.protein and resi 594 and name OD2
set label_position, (-1.6, 1.0, 6), polymer.protein and resi 483 and name CE
set label_position, (1.2, -1.2, 6), polymer.protein and resi 514 and name CD2
orient resn LIG or bolsillo
turn y, -35
turn x, 15
zoom resn LIG or bolsillo, 3
ray 1600, 1250
png ../resultados_prueba/fig_bolsillo.png, dpi=300
