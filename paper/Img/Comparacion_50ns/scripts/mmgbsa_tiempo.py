#!/usr/bin/env python
"""MM/GBSA de 50 ns resuelto en el tiempo: energía por fotograma, medias por segmento
(como el artículo: media ± DE entre segmentos) y descomposición por residuo por bloques."""
import argparse, io, os, json
import numpy as np, pandas as pd

OFFSET = 447
ap = argparse.ArgumentParser()
ap.add_argument("--mmgbsa", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--dt-ns", type=float, default=0.25)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

txt = open(os.path.join(a.mmgbsa, "energias_por_frame.csv")).read().splitlines()
i0 = txt.index("DELTA Energy Terms") + 1
filas = []
for l in txt[i0:]:
    if not l.strip():
        break
    filas.append(l)
E = pd.read_csv(io.StringIO("\n".join(filas)))
E = E.rename(columns={"VDWAALS": "vdw", "EEL": "eel", "EGB": "egb", "ESURF": "esurf", "DELTA TOTAL": "dG",
                      "DELTA G gas": "gas", "DELTA G solv": "solv"})
E["t_ns"] = 0.01 + a.dt_ns * np.arange(len(E))
E = E[["t_ns", "dG", "vdw", "eel", "egb", "esurf", "gas", "solv"]]
E.to_csv(os.path.join(a.out, "mmgbsa_por_frame.csv"), index=False)

R = {"n": len(E)}
for c in ["dG", "vdw", "eel", "egb", "esurf"]:
    R[c] = [round(E[c].mean(), 2), round(E[c].std(), 2)]
# segmentos de 5 ns (20 fotogramas cada uno), como los 20 x 50 ns del artículo
E["seg5"] = np.minimum((E.t_ns // 5).astype(int), 9)
seg5 = E.groupby("seg5")[["dG", "vdw", "eel", "egb", "esurf"]].mean()
R["seg5_dG"] = [round(x, 2) for x in seg5.dG]
R["dG_segmentos_media_sd"] = [round(seg5.dG.mean(), 2), round(seg5.dG.std(), 2)]
E["seg10"] = np.minimum((E.t_ns // 10).astype(int), 4)
seg10 = E.groupby("seg10")[["dG", "vdw", "eel", "egb", "esurf"]].agg(["mean", "std"])
seg10.to_csv(os.path.join(a.out, "mmgbsa_bloques10ns.csv"))
R["seg10_dG"] = [[round(m, 2), round(s, 2)] for m, s in zip(seg10[("dG", "mean")], seg10[("dG", "std")])]
seg5.to_csv(os.path.join(a.out, "mmgbsa_segmentos5ns.csv"))

# descomposición por residuo, por fotograma (sección DELTAS)
lines = open(os.path.join(a.mmgbsa, "decomp_por_frame.csv")).read().splitlines()
j = lines.index("DELTAS:")
k = [n for n in range(j, len(lines)) if lines[n].startswith("DELTA,Total Energy Decomposition")][0]
D = pd.read_csv(io.StringIO("\n".join(lines[k + 1:])))
D.columns = ["frame", "residue", "loc", "internal", "vdw", "elec", "polar", "nopolar", "total"]
D["num"] = D.residue.str.split().str[-1].astype(int)
D["aa"] = D.residue.str.split().str[0]
D["resid"] = np.where(D.aa == "LIG", 0, D.num + OFFSET)
D["t_ns"] = 0.01 + a.dt_ns * (D.frame - 1)
D["bloque10"] = np.minimum((D.t_ns // 10).astype(int), 4)
P = D[D.aa != "LIG"]
media = P.groupby("resid").total.mean().sort_values()
top = media.head(14).index.tolist()
nom = {"HID": "His", "HIE": "His", "HIP": "His", "GLH": "Glu", "ASH": "Asp"}
lab = {r: f"{nom.get(a_, a_.capitalize())}{r}" for r, a_ in zip(P.resid, P.aa)}
tb = P[P.resid.isin(top)].pivot_table(index="resid", columns="bloque10", values="total", aggfunc="mean")
tb["total"] = media[top]
tb = tb.loc[top]; tb.index = [lab[r] for r in tb.index]
tb.to_csv(os.path.join(a.out, "decomp_top_por_bloque10ns.csv"))
comp = P[P.resid.isin(top)].groupby("resid")[["vdw", "elec", "polar", "nopolar", "total"]].mean().loc[top]
comp["sd"] = P[P.resid.isin(top)].groupby("resid").total.std().loc[top]
comp.index = [lab[r] for r in comp.index]
comp.to_csv(os.path.join(a.out, "decomp_top_componentes.csv"))
R["dominantes_le_-1.5"] = [(n, round(v, 2)) for n, v in comp.total.items() if v <= -1.5]
json.dump(R, open(os.path.join(a.out, "mmgbsa_tiempo.json"), "w"), indent=1)
print(json.dumps(R, indent=None))
print(tb.round(2).to_string())
print(comp.round(2).to_string())
