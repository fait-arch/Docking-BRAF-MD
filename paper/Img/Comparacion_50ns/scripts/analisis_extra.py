#!/usr/bin/env python
"""Validación y análisis resuelto en el tiempo de una corrida de 50 ns (producción).

Lee la topología seca y los 10 DCD de producción, comprueba que estén completos y calcula,
fotograma a fotograma: puentes de H ligando-proteína (mismos criterios que 04_analisis_reporte.py:
D-A < 3.5 Å, ángulo D-H-A > 130°), contactos por residuo (< 4 Å entre átomos pesados), distancias
clave, estructura secundaria (DSSP) y enterramiento del ligando. Escribe CSV en --out.
"""
import argparse, glob, json, os
import numpy as np, pandas as pd
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

OFFSET = 447
ap = argparse.ArgumentParser()
ap.add_argument("--tray", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

top = os.path.join(a.tray, "sistema", "complejo_seco.prmtop")
dcds = sorted(glob.glob(os.path.join(a.tray, "md_50ns", "04_prod_*.dcd")))
u = mda.Universe(top, dcds)
V = {"n_dcd": len(dcds), "n_atomos": u.atoms.n_atoms, "n_frames": len(u.trajectory),
     "frames_por_dcd": [len(mda.Universe(top, d).trajectory) for d in dcds]}
nf = len(u.trajectory)
t_ns = 0.01 + 0.01 * np.arange(nf)          # primer fotograma a 10 ps, uno cada 10 ps

lig = u.select_atoms("resname LIG")
ligH = u.select_atoms("resname LIG and not element H")
prot = u.select_atoms("protein and not element H")
resid_prot = prot.resids + OFFSET
res_unicos = np.unique(resid_prot)
resname = {r + OFFSET: n for r, n in zip(prot.residues.resids, prot.residues.resnames)}

# átomos de interés
def sel1(s):
    g = u.select_atoms(s)
    return g
N1 = sel1("resname LIG and name N1")
carbox = u.select_atoms("protein and ((resname ASP ASH and name OD1 OD2) or (resname GLU GLH and name OE1 OE2))")
aromaticos = {}
for rid, rn in [(583, "PHE"), (595, "PHE"), (531, "TRP")]:
    g = u.select_atoms(f"resid {rid-OFFSET} and name CG CD1 CD2 CE1 CE2 CZ" if rn == "PHE" else
                       f"resid {rid-OFFSET} and name CD2 CE2 CE3 CZ2 CZ3 CH2")
    aromaticos[f"{rn.capitalize()}{rid}"] = g
pares = {
    "O1_T529OG1": ("resname LIG and name O1", f"resid {529-OFFSET} and name OG1"),
    "N1_K483NZ": ("resname LIG and name N1", f"resid {483-OFFSET} and name NZ"),
    "O2_K483NZ": ("resname LIG and name O2", f"resid {483-OFFSET} and name NZ"),
    "O2_K591O": ("resname LIG and name O2", f"resid {591-OFFSET} and name O"),
    "N1_D594": ("resname LIG and name N1", f"resid {594-OFFSET} and name OD1 OD2"),
    "N1_E501": ("resname LIG and name N1", f"resid {501-OFFSET} and name OE1 OE2"),
    "lig_T529OG1": ("resname LIG and not element H", f"resid {529-OFFSET} and name OG1"),
    "lig_C532N": ("resname LIG and not element H", f"resid {532-OFFSET} and name N"),
}
pares = {k: (u.select_atoms(s1), u.select_atoms(s2)) for k, (s1, s2) in pares.items()}
pares = {k: v for k, v in pares.items() if len(v[0]) and len(v[1])}

cont = np.zeros((nf, len(res_unicos)), dtype=bool)
dist = {k: np.zeros(nf) for k in pares}
dmin_prot = np.zeros(nf)
carb_min = np.full(nf, np.nan); carb_quien = [""] * nf
cat_pi = {k: np.full(nf, np.nan) for k in aromaticos}
finitos = True
idx_res = {r: i for i, r in enumerate(res_unicos)}
col_de_atomo = np.array([idx_res[r] for r in resid_prot])
for i, ts in enumerate(u.trajectory):
    if not np.all(np.isfinite(u.atoms.positions)):
        finitos = False
    d = distance_array(ligH.positions, prot.positions)
    m = d.min(0)
    dmin_prot[i] = m.min()
    cerca = m < 4.0
    cont[i, np.unique(col_de_atomo[cerca])] = True
    for k, (g1, g2) in pares.items():
        dist[k][i] = distance_array(g1.positions, g2.positions).min()
    if len(N1):
        dc = distance_array(N1.positions, carbox.positions)[0]
        j = int(dc.argmin()); carb_min[i] = dc[j]
        at = carbox[j]; carb_quien[i] = f"{at.resname}{at.resid+OFFSET}"
        for k, g in aromaticos.items():
            cat_pi[k][i] = np.linalg.norm(N1.positions[0] - g.positions.mean(0))
V["coordenadas_finitas"] = finitos
V["dmin_lig_prot_max_A"] = round(float(dmin_prot.max()), 2)

S = pd.DataFrame({"t_ns": t_ns, "dmin_prot": dmin_prot, **{f"d_{k}": v for k, v in dist.items()}})
if len(N1):
    S["d_N1_carbox"] = carb_min; S["carbox_mas_cercano"] = carb_quien
    for k, v in cat_pi.items():
        S[f"d_N1_anillo_{k}"] = v

# puentes de hidrógeno (mismos criterios que el análisis principal)
h = HBA(u, between=["resname LIG", "protein"], d_a_cutoff=3.5, d_h_a_angle_cutoff=130)
h.hydrogens_sel = h.guess_hydrogens("resname LIG or protein")
h.acceptors_sel = h.guess_acceptors("resname LIG or protein")
h.run()
hb = h.results.hbonds   # frame, donor, hydrogen, acceptor, distance(D-A), angle
def etiqueta(ix):
    at = u.atoms[int(ix)]
    return f"{at.resname}{at.resid + (0 if at.resname == 'LIG' else OFFSET)}:{at.name}"
rows = []
for r in hb:
    rows.append({"frame": int(r[0]), "donador": etiqueta(r[1]), "aceptor": etiqueta(r[3]),
                 "d_DA": float(r[4]), "angulo": float(r[5]), "h_idx": int(r[2]), "a_idx": int(r[3])})
HB = pd.DataFrame(rows)
# distancia H···A, como la reporta el artículo
if len(HB):
    dHA = []
    for fr, g in HB.groupby("frame"):
        u.trajectory[fr]
        for _, rr in g.iterrows():
            dHA.append((rr.name, float(np.linalg.norm(u.atoms[rr.h_idx].position - u.atoms[rr.a_idx].position))))
    dHA = dict(dHA); HB["d_HA"] = [dHA[i] for i in HB.index]
    HB["t_ns"] = t_ns[HB.frame.values]
HB.to_csv(os.path.join(a.out, "puentes_H_por_frame.csv"), index=False)
nhb = np.zeros(nf, dtype=int)
if len(HB):
    c = HB.groupby("frame").size(); nhb[c.index.values] = c.values
S["n_puentesH"] = nhb

# bloques de 5 ns: ocupación de contactos y de puentes de H
bloques = np.minimum((t_ns // 5).astype(int), 9)
C = pd.DataFrame(cont, columns=[f"{resname[r].capitalize() if resname[r] not in ('HID','HIE','HIP','GLH','ASH') else {'HID':'His','HIE':'His','HIP':'His','GLH':'Glu','ASH':'Asp'}[resname[r]]}{r}" for r in res_unicos])
C["bloque"] = bloques
occ_bloq = C.groupby("bloque").mean().T
occ_bloq["total"] = C.drop(columns="bloque").mean()
occ_bloq = occ_bloq[occ_bloq.total >= 0.10].sort_values("total", ascending=False)
occ_bloq.to_csv(os.path.join(a.out, "contactos_por_bloque.csv"))
C.drop(columns="bloque").astype(int).to_csv(os.path.join(a.out, "contactos_por_frame.csv.gz"), index=False)

if len(HB):
    HB["bloque"] = bloques[HB.frame.values]
    par = HB.groupby(["donador", "aceptor"])
    tab = []
    for (d_, a_), g in par:
        fr = g.frame.unique()
        fila = {"donador": d_, "aceptor": a_, "ocupacion": len(fr) / nf,
                "d_HA_media": g.d_HA.mean(), "d_HA_sd": g.d_HA.std(), "d_DA_media": g.d_DA.mean()}
        for b in range(10):
            fila[f"b{b}"] = len(set(fr) & set(np.where(bloques == b)[0])) / (bloques == b).sum()
        tab.append(fila)
    pd.DataFrame(tab).sort_values("ocupacion", ascending=False).to_csv(
        os.path.join(a.out, "puentes_H_por_bloque.csv"), index=False)

# estructura secundaria y enterramiento cada 100 ps (mdtraj)
import mdtraj as md
tj = md.load(dcds, top=md.load_prmtop(top), stride=10)
tj = tj.atom_slice(tj.topology.select("protein or resname LIG"))
ss = md.compute_dssp(tj.atom_slice(tj.topology.select("protein")), simplified=True)
li = tj.topology.select("resname LIG")
sasa_c = md.shrake_rupley(tj, mode="atom")[:, li].sum(1)
sasa_l = md.shrake_rupley(tj.atom_slice(li), mode="atom").sum(1)
D = pd.DataFrame({"t_ns": t_ns[::10][:len(ss)], "H": (ss == "H").mean(1) * 100, "E": (ss == "E").mean(1) * 100,
                  "C": (ss == "C").mean(1) * 100, "enterrado_pct": (1 - sasa_c / sasa_l) * 100,
                  "sasa_lig_A2": sasa_c * 100})
D.to_csv(os.path.join(a.out, "dssp_sasa_cada100ps.csv"), index=False)
np.save(os.path.join(a.out, "dssp.npy"), ss)
S.to_csv(os.path.join(a.out, "series_extra.csv"), index=False)
json.dump(V, open(os.path.join(a.out, "validacion.json"), "w"), indent=1)
print(json.dumps(V))
