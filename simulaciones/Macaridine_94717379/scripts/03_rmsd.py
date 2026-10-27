#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
03_rmsd.py — RMSD de la trayectoria (proteína Cα y ligando) y gráfica.

  * Referencia: estructura minimizada (01_min.pdb), es decir, el complejo de docking
    relajado antes de empezar la dinámica.
  * Proteína: RMSD de Cα tras superponer sobre Cα.
  * Ligando: RMSD de átomos pesados del ligando SIN volver a superponer (sólo el
    ajuste sobre Cα de la proteína), es decir, cuánto se desplaza dentro del bolsillo.

Uso:
  python 03_rmsd.py --sistema ../sistema --md ../md_prueba --out ../resultados
"""
import argparse, os, warnings
import numpy as np
warnings.filterwarnings("ignore")

FASES = [("02_cal", "Calentamiento"), ("03_eq", "Equilibrado"), ("04_prod_*", "Producción")]

# Paleta de referencia (skill dataviz): superficie, tinta y dos series categóricas
SUP, TXT, TXT2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
AZUL, NARANJA = "#2a78d6", "#eb6834"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sistema", default="sistema")
    ap.add_argument("--md", default="md")
    ap.add_argument("--out", default="resultados")
    ap.add_argument("--titulo", default="Macaridine (94717379) – BRAF V600E")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    import pandas as pd

    top = os.path.join(a.sistema, "sistema.prmtop")
    ref = mda.Universe(top, os.path.join(a.md, "01_min.pdb"), format="PDB")

    filas = [{"fase": "Minimizada", "tiempo_ps": 0.0, "rmsd_CA_A": 0.0, "rmsd_ligando_A": 0.0}]
    limites, t0 = [], 0.0
    import glob
    for pref, nombre in FASES:
        dcds = sorted(glob.glob(os.path.join(a.md, pref + ".dcd")))
        dcds = [d for d in dcds if os.path.getsize(d) > 0]
        if not dcds:
            continue
        # la trayectoria puede contener el sistema completo o sólo proteína+ligando
        from MDAnalysis.lib.formats.libdcd import DCDFile
        with DCDFile(dcds[0]) as fh:
            n_dcd = fh.header["natoms"]
        top_i = top if n_dcd == len(ref.atoms) else os.path.join(a.sistema, "complejo_seco.prmtop")
        u = mda.Universe(top_i, dcds)
        dt = u.trajectory.dt  # ps entre fotogramas
        ref_i = ref.atoms if top_i == top else ref.atoms[:n_dcd]
        R = rms.RMSD(u, ref_i, select="protein and name CA",
                     groupselections=["resname LIG and not element H"]).run()
        res = R.results.rmsd
        t = t0 + (np.arange(len(res)) + 1) * dt
        for ti, row in zip(t, res):
            filas.append({"fase": nombre, "tiempo_ps": round(float(ti), 3),
                          "rmsd_CA_A": float(row[2]), "rmsd_ligando_A": float(row[3])})
        limites.append((t0, float(t[-1]), nombre))
        t0 = float(t[-1])

    df = pd.DataFrame(filas)
    csv = os.path.join(a.out, "rmsd.csv")
    df.to_csv(csv, index=False)
    prod = df[df.fase == "Producción"]
    resumen = df[df.fase != "Minimizada"].groupby("fase", sort=False)[["rmsd_CA_A", "rmsd_ligando_A"]].agg(["mean", "std", "max"]).round(2)
    print(resumen.to_string())
    resumen.to_csv(os.path.join(a.out, "rmsd_resumen.csv"))

    # ------------------------------------------------------------------ gráfica
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.edgecolor": GRID, "axes.labelcolor": TXT2,
                         "xtick.color": TXT2, "ytick.color": TXT2,
                         "text.color": TXT})
    # eje de tiempo en ns si la trayectoria es larga
    esc, unidad = (1000.0, "ns") if t0 > 5000 else (1.0, "ps")
    df["t_plot"] = df.tiempo_ps / esc
    limites = [(i / esc, f / esc, n) for i, f, n in limites]
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 6.2), sharex=True, dpi=200,
                             gridspec_kw={"hspace": 0.32})
    fig.patch.set_facecolor(SUP)
    series = [("rmsd_CA_A", "Proteína — Cα", AZUL), ("rmsd_ligando_A", "Ligando — átomos pesados (ajuste sobre Cα)", NARANJA)]
    for ax, (col, etiqueta, color) in zip(axes, series):
        ax.set_facecolor(SUP)
        for i, (ini, fin, nombre) in enumerate(limites):
            if nombre != "Producción":
                ax.axvspan(ini, fin, color="#f0efec", lw=0, zorder=0)
            if i > 0:
                ax.axvline(ini, color=GRID, lw=1, zorder=1)
        ax.plot(df.t_plot, df[col], color=color, lw=1.6, solid_capstyle="round", zorder=3)
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_ylabel("RMSD (Å)")
        ax.set_ylim(0, max(df[col].max() * 1.25, 0.5))
        ax.set_title(etiqueta, loc="left", fontsize=10.5, color=TXT, fontweight="bold", pad=6)
        if len(prod):
            m = prod[col].mean()
            ax.set_title(f"media en producción: {m:.2f} ± {prod[col].std():.2f} Å", loc="right",
                         fontsize=9, color=TXT2, pad=6)
    ymax = axes[0].get_ylim()[1]
    for ini, fin, nombre in limites:
        axes[0].text((ini + fin) / 2, ymax * 0.97, nombre, ha="center", va="top", fontsize=8.5, color=TXT2)
    axes[-1].set_xlabel(f"Tiempo ({unidad})")
    axes[-1].set_xlim(0, df.t_plot.max())
    fig.suptitle(a.titulo + " · RMSD respecto a la estructura minimizada",
                 x=0.075, ha="left", fontsize=12, fontweight="bold", color=TXT, y=0.985)
    png = os.path.join(a.out, "rmsd.png")
    fig.savefig(png, bbox_inches="tight", facecolor=SUP)
    print(f"\nGráfica: {png}\nDatos:   {csv}")


if __name__ == "__main__":
    main()
