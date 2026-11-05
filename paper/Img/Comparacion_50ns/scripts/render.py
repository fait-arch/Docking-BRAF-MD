# Paneles 3D del sitio de unión para las conformaciones de 50 ns.  Uso: pymol -cq render.py
from pymol import cmd
import os
os.chdir("/home/claude/report/an50/pymol")

cmd.set("ray_opaque_background", 1); cmd.bg_color("white")
cmd.set("ray_shadows", 0); cmd.set("antialias", 2); cmd.set("ray_trace_mode", 1)
cmd.set("ray_trace_color", "0x444444"); cmd.set("ray_trace_gain", 0.12)
cmd.set("cartoon_fancy_helices", 1); cmd.set("cartoon_side_chain_helper", 1); cmd.set("valence", 0)
cmd.set("cartoon_transparency", 0.82); cmd.set("stick_radius", 0.17)
cmd.set("label_size", 17); cmd.set("float_labels", 1); cmd.set("label_color", "black")
cmd.set("label_font_id", 7); cmd.set("ray_label_specular", 0)
cmd.set("dash_color", "black"); cmd.set("dash_radius", 0.06); cmd.set("dash_gap", 0.25)
cmd.set("dash_width", 2.5); cmd.set("label_distance_digits", 1)
cmd.set_color("c_gris", [0.80, 0.80, 0.78]); cmd.set_color("c_res", [0.62, 0.66, 0.70])
cmd.set_color("c_azul", [0.165, 0.471, 0.839]); cmd.set_color("c_naranja", [0.922, 0.408, 0.204])
cmd.set_color("c_amar", [0.929, 0.631, 0.0]); cmd.set_color("c_mag", [0.910, 0.482, 0.643])

RES = {"M": [463, 471, 481, 483, 514, 529, 531, 532, 583, 591, 594],
       "P": [463, 514, 529, 532, 536, 581, 583, 591, 594]}
LAB = {"M": [471, 483, 514, 529, 532, 583, 591, 594],
       "P": [463, 514, 529, 532, 536, 581, 583, 591, 594]}
ETQ = {463: "CD1", 471: "CG1", 481: "CB", 483: "NZ", 514: "CD1", 529: "OG1", 531: "CZ2", 532: "C",
       535: "OG", 536: "OG", 580: "ND2", 581: "ND2", 583: "CZ", 591: "CD", 594: "CG"}
NOM = {463: "Ile463", 471: "Val471", 481: "Ala481", 483: "Lys483", 514: "Leu514", 529: "Thr529", 531: "Trp531",
       532: "Cys532", 535: "Ser535", 536: "Ser536", 580: "Asn580", 581: "Asn581", 583: "Phe583", 591: "Lys591",
       594: "Asp594"}
PARES = {"M": [("LIG", "O1", 529, "OG1"), ("LIG", "N1", 483, "NZ"), ("LIG", "O1", 483, "NZ")],
         "P": [("LIG", "N1", 594, "OD1"), ("LIG", "N1", 594, "OD2"), ("LIG", "O3", 463, "O")]}

for f in ["M_C", "M_B", "M_A", "P_C", "P_A", "P_F"]:
    cmd.load(f + ".pdb", f)
cmd.remove("hydro and not (resn LIG and neighbor (elem N+O))")
# vista común: bolsillo de la conformación inicial de Macaridine
cmd.hide("everything")
cmd.orient("M_C and resn LIG")
cmd.zoom("(M_C and resn LIG) or (M_C and resi 483+529+532+583+591+594+463 and not name N+C+O) or (P_F and resn LIG)", 1.5)
cmd.turn("x", -25)
cmd.clip("slab", 38)
vista = cmd.get_view()


def panel(obj, sis):
    cmd.hide("everything"); cmd.delete("d_*")
    cmd.set_view(vista)
    cmd.show("cartoon", f"{obj} and polymer.protein"); cmd.color("c_gris", f"{obj} and polymer.protein")
    cmd.color("c_amar", f"{obj} and polymer.protein and resi 464-471")
    cmd.color("c_mag", f"{obj} and polymer.protein and resi 594-623")
    sel = f"{obj} and resi {'+'.join(map(str, RES[sis]))} and not name N+C+O"
    cmd.show("sticks", sel)
    cmd.show("sticks", f"{obj} and resi 532+593+591 and name N+CA+C+O")
    cmd.color("c_res", f"({sel} or ({obj} and resi 532+593+591)) and elem C")
    cmd.color("red", f"{obj} and polymer.protein and elem O"); cmd.color("blue", f"{obj} and polymer.protein and elem N")
    cmd.color("yellow", f"{obj} and polymer.protein and elem S")
    cmd.show("sticks", f"{obj} and resn LIG")
    cmd.set("stick_radius", 0.26, f"{obj} and resn LIG")
    cmd.color("c_naranja" if sis == "P" else "c_azul", f"{obj} and resn LIG and elem C")
    cmd.color("red", f"{obj} and resn LIG and elem O"); cmd.color("blue", f"{obj} and resn LIG and elem N")
    cmd.color("white", f"{obj} and resn LIG and elem H")
    cmd.hide("cartoon", f"{obj} and not polymer.protein")
    for i, (_, la, r, ra) in enumerate(PARES[sis]):
        a1, a2 = f"{obj} and resn LIG and name {la}", f"{obj} and resi {r} and name {ra}"
        if cmd.count_atoms(a1) and cmd.count_atoms(a2) and cmd.get_distance(a1, a2) <= 3.6:
            cmd.distance(f"d_{i}", a1, a2); cmd.color("black", f"d_{i}")
    cmd.hide("labels")
    for r in LAB[sis]:
        cmd.label(f"{obj} and resi {r} and name {ETQ[r]}", f'"{NOM[r]}"')
    cmd.set("label_position", (0, 0, 8)); cmd.set("label_color", "black", obj)
    cmd.ray(1100, 900)
    cmd.png(f"panel_{obj}.png", dpi=300)


for obj in ["M_C", "M_B", "M_A"]:
    panel(obj, "M")
for obj in ["P_C", "P_A", "P_F"]:
    panel(obj, "P")
