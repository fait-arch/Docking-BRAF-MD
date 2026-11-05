#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
05_mmgbsa_figura.py — resume la salida de MMPBSA.py (MM/GBSA + descomposición por residuo)
y genera fig_mmgbsa.pdf/.png y mmgbsa_resumen.json.

  python 05_mmgbsa_figura.py --mmgbsa ../analisis/mmgbsa --fig <carpeta_figuras>
"""
import argparse, csv, io, json, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OFFSET = 447
TXT, TXT2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
C1, C2, C2L = "#2a78d6", "#eb6834", "#f2b199"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#c9c8c2",
                     "axes.labelcolor": TXT2, "xtick.color": TXT2, "ytick.color": TXT2,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                     "savefig.bbox": "tight", "pdf.fonttype": 42})
NOMBRE = {"HID": "His", "HIE": "His", "HIP": "His", "GLH": "Glu", "ASH": "Asp", "CYX": "Cys", "LYN": "Lys"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mmgbsa", required=True)
    ap.add_argument("--fig", required=True)
    ap.add_argument("--dt-ps", type=float, default=1.0, help="ps entre los fotogramas analizados")
    a = ap.parse_args()
    d = a.mmgbsa

    # --- resumen global
    txt = open(os.path.join(d, "FINAL_MMGBSA.dat")).read()
    bloque = txt[txt.index("Differences (Complex - Receptor - Ligand)"):]
    val = {}
    for k in ("VDWAALS", "EEL", "EGB", "ESURF", "DELTA G gas", "DELTA G solv", "DELTA TOTAL"):
        m = re.search(rf"^{re.escape(k)}\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", bloque, re.M)
        val[k] = [float(m.group(1)), float(m.group(2)), float(m.group(3))]

    # --- ΔG por fotograma
    lineas = open(os.path.join(d, "energias_por_frame.csv")).read().splitlines()
    i0 = lineas.index("DELTA Energy Terms") + 1
    filas = []
    for l in lineas[i0:]:
        if not l.strip() or not l[0].isdigit() and not l.startswith("Frame"):
            break
        filas.append(l)
    dfr = pd.read_csv(io.StringIO("\n".join(filas)))
    dg = dfr["DELTA TOTAL"].values

    # --- descomposición por residuo (bloque DELTAS, TOTAL)
    dl = open(os.path.join(d, "FINAL_DECOMP.dat")).read().splitlines()
    j = dl.index("DELTAS:")
    rows = []
    for l in dl[j + 4:]:
        if not l.strip():
            break
        p = next(csv.reader([l]))
        if len(p) < 18:
            break
        resn, num = p[0].split()
        rows.append({"res": resn, "num": int(num), "vdw": float(p[5]), "elec": float(p[8]),
                     "polar": float(p[11]), "nopolar": float(p[14]), "total": float(p[17]), "sd": float(p[18])})
    dec = pd.DataFrame(rows)
    dec = dec[dec.res != "LIG"].copy()
    dec["residuo"] = dec.num + OFFSET
    dec["etiqueta"] = [f"{NOMBRE.get(r, r.capitalize())}{n}" for r, n in zip(dec.res, dec.residuo)]
    dec.sort_values("total").to_csv(os.path.join(d, "descomposicion_por_residuo.csv"), index=False)
    top = dec.sort_values("total").head(10)

    res = {"dG": val["DELTA TOTAL"], "vdw": val["VDWAALS"], "eel": val["EEL"], "egb": val["EGB"],
           "esurf": val["ESURF"], "gas": val["DELTA G gas"], "solv": val["DELTA G solv"],
           "n_frames": int(len(dg)),
           "top": [(e, round(t, 2), round(v, 2), round(el, 2)) for e, t, v, el in
                   zip(top.etiqueta, top.total, top.vdw, top.elec)],
           "dominantes": [e for e, t in zip(top.etiqueta, top.total) if t <= -1.5]}
    json.dump(res, open(os.path.join(d, "mmgbsa_resumen.json"), "w"), indent=2, ensure_ascii=False)
    print(json.dumps(res, indent=2, ensure_ascii=False))

    # --- figura
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"width_ratios": [1.15, 1], "wspace": 0.38})
    t = np.arange(1, len(dg) + 1) * a.dt_ps
    esc, uni = (1000.0, "ns") if t[-1] > 5000 else (1.0, "ps")
    t = t / esc
    ax1.plot(t, dg, color=C1, lw=1.2)
    ax1.axhline(dg.mean(), color=TXT2, lw=0.9, ls=(0, (3, 3)))
    ax1.text(0.98, 0.03, f"media {dg.mean():.1f} ± {dg.std():.1f} kcal/mol".replace("-", "\u2212"),
             transform=ax1.transAxes, ha="right", va="bottom", fontsize=7.5, color=TXT2,
             bbox=dict(fc="white", ec="none", pad=1.5))
    ax1.set_xlabel(f"Tiempo de producción ({uni})"); ax1.set_ylabel(r"$\Delta G_{\mathrm{unión}}$ (kcal/mol)")
    ax1.set_xlim(0, t[-1] * 1.01)
    y = np.arange(len(top))[::-1]
    col = [C2 if v <= -1.5 else C2L for v in top.total]
    ax2.barh(y, top.total, xerr=top.sd, color=col, height=0.68, error_kw={"elinewidth": 0.7, "ecolor": TXT2, "capsize": 0})
    ax2.set_yticks(y); ax2.set_yticklabels(top.etiqueta, fontsize=8)
    ax2.axvline(-1.5, color=TXT2, lw=0.8, ls=(0, (3, 3)))
    ax2.axvline(0, color="#c9c8c2", lw=0.8)
    ax2.set_xlabel("Contribución (kcal/mol)")
    ax2.grid(axis="y", visible=False)
    for s in ("pdf", "png"):
        fig.savefig(os.path.join(a.fig, f"fig_mmgbsa.{s}"), dpi=170, facecolor="white")


if __name__ == "__main__":
    main()
