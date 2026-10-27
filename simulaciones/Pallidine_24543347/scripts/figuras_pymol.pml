# Figuras 3D del reporte de Pallidine (PyMOL). Ejecutar desde scripts/:  pymol -cq figuras_pymol.pml
# Estructura final de la producción (cx) y pose de partida minimizada (ini) superpuestas sobre Cα.
load ../md_prueba/04_prod_004.pdb, cx
load ../md_prueba/01_min.pdb, ini
remove not (polymer.protein or resn LIG)
alter polymer.protein, resi=str(int(resi)+447)
sort
super ini and polymer.protein and name CA, cx and polymer.protein and name CA
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
set_color c_gris_osc, [0.55, 0.55, 0.53]
hide everything
show cartoon, cx and polymer.protein
color c_gris, cx and polymer.protein
color c_amar, cx and polymer.protein and resi 464-471
color c_viol, cx and polymer.protein and resi 492-507
color c_aqua, cx and polymer.protein and resi 529-535
color c_mag, cx and polymer.protein and resi 594-623
show spheres, cx and resn LIG and not hydro
color c_naranja, cx and resn LIG and elem C
color red, cx and resn LIG and elem O
color blue, cx and resn LIG and elem N
orient cx and polymer.protein
turn z, -90
zoom cx and polymer.protein, -3
ray 1500, 1500
png ../resultados_prueba/fig_complejo_general.png, dpi=300

# ------------------------------------------------ acercamiento al bolsillo: dos paneles, misma vista
# A = pose de docking minimizada (ini), B = final de la producción (cx)
hide spheres
set cartoon_color, c_gris
set cartoon_transparency, 0.86
set stick_radius, 0.17
set label_size, -1.5
set float_labels, 1
set ray_label_specular, 0
set label_color, 0x111111
set label_outline_color, white
set label_font_id, 7
set label_position, (0, 0, 6)
set dash_color, 0x111111
set dash_radius, 0.07
set dash_gap, 0.22
set_color c_lig_ini, [0.55, 0.55, 0.53]
select resB, polymer.protein and resi 463+483+514+529+530+532+536+580+581+583+591+594 and not name N+C+O
select resB_bb, polymer.protein and resi 532+593 and name N+CA+C+O
orient cx and resn LIG
zoom (cx and resn LIG) or (cx and resB) or (ini and resn LIG), 2
clip slab, 40
view v1, store

python
etq = [(483,"NZ","Lys483",(0.6,-1.2,6)), (514,"CD1","Leu514",(0.9,-1.5,15)), (529,"OG1","Thr529",(-3.2,0.6,15)),
       (530,"NE2","Gln530",(1.6,0.9,6)), (532,"SG","Cys532",(1.4,0.8,6)), (536,"OG","Ser536",(0,-0.9,6)),
       (580,"ND2","Asn580",(0,-0.9,6)), (581,"ND2","Asn581",(-1.4,0,6)), (583,"CZ","Phe583",(2.4,0.8,6)),
       (591,"NZ","Lys591",(-1.8,0.9,6)), (594,"CG","Asp594",(-1.8,0,6)), (463,"CD1","Ile463",(1.2,0,6))]
def panel(obj, ligcolor):
    cmd.hide("sticks"); cmd.hide("labels"); cmd.delete("hb*")
    cmd.show("sticks", f"{obj} and (resB or resB_bb) and not hydro")
    cmd.color("c_azul", f"{obj} and (resB or resB_bb) and elem C")
    cmd.color("red", f"{obj} and (resB or resB_bb) and elem O")
    cmd.color("blue", f"{obj} and (resB or resB_bb) and elem N")
    cmd.color("yellow", f"{obj} and resB and elem S")
    cmd.show("sticks", f"{obj} and resn LIG and not hydro")
    cmd.show("sticks", f"{obj} and resn LIG and hydro and neighbor elem O")
    cmd.color(ligcolor, f"{obj} and resn LIG and elem C")
    cmd.color("red", f"{obj} and resn LIG and elem O")
    cmd.color("blue", f"{obj} and resn LIG and elem N")
    cmd.color("white", f"{obj} and resn LIG and hydro")
    cmd.set("stick_radius", 0.27, f"{obj} and resn LIG")
    for r, at, txt, pos in etq:
        s = f"{obj} and polymer.protein and resi {r} and name {at}"
        cmd.label(s, repr(txt)); cmd.set("label_position", pos, s)
    cmd.view("v1", "recall")
python end

python
panel("ini", "c_lig_ini")
cmd.distance("hb1", "ini and resn LIG and name O1", "ini and polymer.protein and resi 532 and name N", 3.7)
cmd.hide("labels", "hb1")
python end
ray 1500, 1250
png ../resultados_prueba/fig_bolsillo_A.png, dpi=300

python
panel("cx", "c_naranja")
python end
ray 1500, 1250
png ../resultados_prueba/fig_bolsillo_B.png, dpi=300
