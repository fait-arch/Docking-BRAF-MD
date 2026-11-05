#!/usr/bin/env python
"""Figuras comparativas de las corridas de 50 ns (Macaridine y Pallidine NH+) para simulacion.tex."""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

BASE = "/home/claude/report/an50"
RES = "/mnt/user-data/uploads/Github/Docking-BRAF-MD/resultados_50ns"
FIG = sys.argv[1] if len(sys.argv) > 1 else os.path.join(BASE, "fig")
os.makedirs(FIG, exist_ok=True)
SIS = {"M": ("Macaridine", "Macaridine_94717379"), "P": ("Pallidine NH$^+$", "Pallidine_24543347_protonada")}

SUP, TXT, TXT2, GRID, BAND = "#ffffff", "#0b0b0b", "#52514e", "#e4e3df", "#f0efec"
C1, C2, C3, C4, C5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
COL = {"M": C1, "P": C2}
CLU = [C3, C5, C4, "#6250d6", "#8c8b86"]
REGIONES = [(464, 471, "Lazo P"), (492, 507, "Hélice αC"), (529, 535, "Bisagra"),
            (574, 581, "Lazo cat."), (594, 596, "DFG"), (597, 623, "Lazo de activación")]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 9,
    "axes.edgecolor": "#c9c8c2", "axes.labelcolor": TXT2, "xtick.color": TXT2, "ytick.color": TXT2,
    "text.color": TXT, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "lines.linewidth": 1.0, "savefig.bbox": "tight", "savefig.pad_inches": 0.04, "pdf.fonttype": 42,
    "legend.frameon": False, "legend.fontsize": 7.8,
})


def guardar(fig, nombre):
    fig.savefig(os.path.join(FIG, nombre + ".pdf"), facecolor=SUP)
    fig.savefig(os.path.join(FIG, nombre + ".png"), dpi=170, facecolor=SUP)
    plt.close(fig)


def suave(y, n=50):
    return pd.Series(y).rolling(n, center=True, min_periods=1).mean().values


def titulo(ax, txt, der=None):
    ax.set_title(txt, loc="left", fontweight="bold")
    if der:
        ax.set_title(der, loc="right", color=TXT2, fontsize=7.8)


D = {}
for k, (nom, car) in SIS.items():
    s = pd.read_csv(os.path.join(RES, car, "analisis_50ns", "series_estructurales.csv"))
    s = s[s.fase == "Producción"].reset_index(drop=True)
    s["t_ns"] = 0.01 + 0.01 * np.arange(len(s))
    x = pd.read_csv(os.path.join(BASE, "out", k, "series_extra.csv"))
    clu = os.path.join(BASE, "cluster", k, "conf_por_frame.csv")
    D[k] = {"s": s, "x": x,
            "rmsf": pd.read_csv(os.path.join(RES, car, "analisis_50ns", "rmsf.csv")),
            "gb": pd.read_csv(os.path.join(BASE, "out", f"{k}_gb", "mmgbsa_por_frame.csv")),
            "gbseg": pd.read_csv(os.path.join(BASE, "out", f"{k}_gb", "mmgbsa_segmentos5ns.csv")),
            "dec": pd.read_csv(os.path.join(BASE, "out", f"{k}_gb", "decomp_top_componentes.csv"), index_col=0),
            "occ": pd.read_csv(os.path.join(BASE, "out", k, "contactos_por_bloque.csv"), index_col=0),
            "ss": pd.read_csv(os.path.join(BASE, "out", k, "dssp_sasa_cada100ps.csv")),
            "clu": pd.read_csv(clu) if os.path.exists(clu) else None}

# ---------------------------------------------------------------- 1) estabilidad (como la Fig. 6 del artículo)
fig, axs = plt.subplots(3, 2, figsize=(7.2, 5.6), sharex=True, gridspec_kw={"hspace": 0.42, "wspace": 0.22})
for j, k in enumerate(["M", "P"]):
    s, c = D[k]["s"], D[k]["clu"]
    ax = axs[0, j]
    if c is not None:
        for i, lab in enumerate(sorted(c.conf.unique())):
            m = (c.conf == lab).values
            ax.scatter(s.t_ns[m], s.rmsd_CA[m], s=0.6, color=CLU[i], lw=0, rasterized=True)
        hand = [Line2D([], [], marker="o", ls="", color=CLU[i], ms=4, label=f"Conf. {lab}")
                for i, lab in enumerate(sorted(c.conf.unique()))]
        ax.legend(handles=hand, loc="lower right", ncol=3, handletextpad=0.1, columnspacing=0.6)
    else:
        ax.plot(s.t_ns, s.rmsd_CA, color=COL[k], lw=0.6)
    ax.set_ylim(0, 3.2); ax.set_ylabel("RMSD Cα (Å)")
    titulo(ax, SIS[k][0], f"{s.rmsd_CA.mean():.2f} ± {s.rmsd_CA.std():.2f} Å")
    ax = axs[1, j]
    ax.plot(s.t_ns, s.rmsd_lig, color=COL[k], lw=0.5, alpha=0.45)
    ax.plot(s.t_ns, suave(s.rmsd_lig), color=COL[k], lw=1.3)
    ax.set_ylim(0, 8); ax.set_ylabel("RMSD ligando (Å)")
    titulo(ax, "Ligando", f"{s.rmsd_lig.mean():.2f} ± {s.rmsd_lig.std():.2f} Å")
    ax = axs[2, j]
    ax.plot(s.t_ns, s.rg, color=TXT2, lw=0.5, alpha=0.5)
    ax.plot(s.t_ns, suave(s.rg), color=TXT, lw=1.2)
    ax.set_ylim(18.8, 19.9); ax.set_ylabel("$R_g$ proteína (Å)")
    ax.set_xlabel("Tiempo de producción (ns)")
    titulo(ax, "Radio de giro", f"{s.rg.mean():.2f} ± {s.rg.std():.2f} Å")
    ax.set_xlim(0, 50)
guardar(fig, "fig50_estabilidad")

# ---------------------------------------------------------------- 2) RMSF (como la Fig. 8 del artículo)
fig, ax = plt.subplots(figsize=(7.2, 2.6))
ymax = 4.0
for i, f, n in REGIONES:
    ax.axvspan(i - 0.5, f + 0.5, color=BAND, lw=0, zorder=0)
    xx, yy = ((i + f) / 2 + 9, ymax * 0.985) if n == "Lazo de activación" else ((i + f) / 2, ymax * 0.985)
    if n == "DFG":
        yy = 2.6
    if n == "Lazo P":
        xx = f + 4
    ax.text(xx, yy, n, ha="center", va="top", fontsize=6.8, color=TXT2, rotation=90 if f - i < 9 else 0)
for k in ["M", "P"]:
    r = D[k]["rmsf"]
    ax.plot(r.residuo, r.rmsf_A, color=COL[k], lw=1.0, label=f"{SIS[k][0]} (media {r.rmsf_A.mean():.2f} Å)")
ax.axhline(1.0, color=TXT2, lw=0.6, ls=(0, (3, 3)))
ax.set_xlim(448, 723); ax.set_ylim(0, ymax)
ax.set_xlabel("Residuo de BRAF"); ax.set_ylabel("RMSF Cα (Å)")
ax.legend(loc="center left", bbox_to_anchor=(0.66, 0.62))
guardar(fig, "fig50_rmsf")

# ---------------------------------------------------------------- 3) puentes de H y distancias clave (Fig. 9)
fig, axs = plt.subplots(2, 2, figsize=(7.2, 4.4), sharex=True, gridspec_kw={"hspace": 0.45, "wspace": 0.25})
xm, xp = D["M"]["x"], D["P"]["x"]
ax = axs[0, 0]
ax.plot(xm.t_ns, xm.d_O1_T529OG1, color=C1, lw=0.4, alpha=0.4)
ax.plot(xm.t_ns, suave(xm.d_O1_T529OG1), color=C1, lw=1.3, label="O1(N–OH)···OG1 Thr529")
ax.plot(xm.t_ns, suave(xm.d_N1_K483NZ), color=C3, lw=1.1, label="N1···NZ Lys483")
ax.axhline(3.5, color=TXT2, lw=0.7, ls=(0, (3, 3)))
ax.set_ylim(2, 9); ax.set_ylabel("Distancia (Å)")
titulo(ax, "Macaridine: anclajes polares")
ax.legend(loc="upper left", ncol=1)
ax = axs[1, 0]
ax.fill_between(xm.t_ns, suave(xm.n_puentesH, 50), color=C1, alpha=0.25, lw=0)
ax.plot(xm.t_ns, suave(xm.n_puentesH, 50), color=C1, lw=1.2)
ax.set_ylim(0, 1.6); ax.set_ylabel("Puentes de H (media móvil)")
ax.set_xlabel("Tiempo de producción (ns)")
titulo(ax, "Puentes de H", f"media {xm.n_puentesH.mean():.2f}")
ax = axs[0, 1]
ax.plot(xp.t_ns, xp.d_N1_carbox, color=C2, lw=0.4, alpha=0.4)
ax.plot(xp.t_ns, suave(xp.d_N1_carbox), color=C2, lw=1.3, label="N$^+$···carboxilato más cercano")
ax.plot(xp.t_ns, suave(xp["d_N1_anillo_Phe583"]), color=C4, lw=1.1, label="N$^+$···anillo Phe583")
ax.axhline(4.0, color=TXT2, lw=0.7, ls=(0, (3, 3)))
ax.set_ylim(2, 9)
titulo(ax, "Pallidine NH$^+$: grupo catiónico")
ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.64))
ax = axs[1, 1]
ax.fill_between(xp.t_ns, suave(xp.n_puentesH, 50), color=C2, alpha=0.25, lw=0)
ax.plot(xp.t_ns, suave(xp.n_puentesH, 50), color=C2, lw=1.2)
ax.set_ylim(0, 1.6); ax.set_xlabel("Tiempo de producción (ns)")
titulo(ax, "Puentes de H", f"media {xp.n_puentesH.mean():.2f}")
for ax in axs.flat:
    ax.set_xlim(0, 50)
guardar(fig, "fig50_puentesH")

# ---------------------------------------------------------------- 4) contactos por bloque de 5 ns
fig, axs = plt.subplots(1, 2, figsize=(7.2, 4.3), gridspec_kw={"wspace": 0.42})
from matplotlib.colors import LinearSegmentedColormap
for ax, k in zip(axs, ["M", "P"]):
    o = D[k]["occ"]; o = o[o.total >= 0.25].head(18)
    cmap = LinearSegmentedColormap.from_list("c", ["#ffffff", COL[k]])
    m = o[[str(b) for b in range(10)]].values * 100
    im = ax.imshow(m, aspect="auto", cmap=cmap, vmin=0, vmax=100)
    ax.set_yticks(range(len(o))); ax.set_yticklabels(o.index, fontsize=7.5)
    ax.set_xticks(range(10)); ax.set_xticklabels([f"{5*b}–{5*b+5}" for b in range(10)], rotation=60, fontsize=6.8)
    ax.set_xlabel("Bloque de producción (ns)"); ax.grid(False)
    for (i, j), v in np.ndenumerate(m):
        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=5.6, color="white" if v > 60 else TXT2)
    titulo(ax, SIS[k][0])
    for sp in ax.spines.values():
        sp.set_visible(False)
cb = fig.colorbar(im, ax=axs, fraction=0.025, pad=0.02); cb.set_label("% de fotogramas en contacto (< 4 Å)")
cb.outline.set_visible(False)
guardar(fig, "fig50_contactos")

# ---------------------------------------------------------------- 5) MM/GBSA (Tabla 4 / Fig. 10 del artículo)
fig = plt.figure(figsize=(7.2, 5.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.5, wspace=0.42)
ax = fig.add_subplot(gs[0, :])
for k in ["M", "P"]:
    g, sg = D[k]["gb"], D[k]["gbseg"]
    ax.plot(g.t_ns, g.dG, color=COL[k], lw=0.6, alpha=0.45)
    xs = np.repeat(np.arange(0, 55, 5), 2)[1:-1]
    ys = np.repeat(sg.dG.values, 2)
    ax.plot(xs, ys, color=COL[k], lw=2.0, label=f"{SIS[k][0]}: {g.dG.mean():.1f} ± {g.dG.std():.1f} kcal/mol")
ax.set_xlim(0, 50); ax.set_xlabel("Tiempo de producción (ns)"); ax.set_ylabel("Δ$G_{\\mathrm{unión}}$ (kcal/mol)")
ax.legend(loc="upper right", ncol=2)
titulo(ax, "Energía libre de unión por fotograma y media por segmento de 5 ns")
for col, k in enumerate(["M", "P"]):
    ax = fig.add_subplot(gs[1, col])
    d = D[k]["dec"].head(10).iloc[::-1]
    yy = np.arange(len(d))
    left = np.zeros(len(d))
    ax.barh(yy, d.total, color=[COL[k] if v <= -1.5 else ("#9cc0ec" if k == "M" else "#f2b199") for v in d.total],
            height=0.68, xerr=d.sd, error_kw={"elinewidth": 0.7, "ecolor": TXT2, "capsize": 0})
    ax.axvline(-1.5, color=TXT2, lw=0.7, ls=(0, (3, 3)))
    ax.set_yticks(yy); ax.set_yticklabels(d.index, fontsize=7.6)
    ax.set_xlim(-5.2, 0.3); ax.set_xlabel("Contribución (kcal/mol)")
    ax.grid(axis="y", visible=False)
    titulo(ax, SIS[k][0])
guardar(fig, "fig50_mmgbsa")

# ---------------------------------------------------------------- 6) estructura secundaria y enterramiento
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.3), gridspec_kw={"wspace": 0.3})
for k in ["M", "P"]:
    ss = D[k]["ss"]
    axs[0].plot(ss.t_ns, suave(ss.H, 10), color=COL[k], lw=1.2, label=f"{SIS[k][0]} hélice")
    axs[0].plot(ss.t_ns, suave(ss.E, 10), color=COL[k], lw=1.2, ls=(0, (3, 2)), label=f"{SIS[k][0]} lámina β")
    axs[1].plot(ss.t_ns, suave(ss.enterrado_pct, 10), color=COL[k], lw=1.2, label=SIS[k][0])
axs[0].set_ylim(0, 50); axs[0].set_ylabel("% de residuos"); titulo(axs[0], "Estructura secundaria (DSSP)")
axs[0].legend(loc="center right", fontsize=6.8, ncol=1)
axs[1].set_ylim(80, 100); axs[1].set_ylabel("% de SASA enterrada"); titulo(axs[1], "Enterramiento del ligando")
axs[1].legend(loc="lower right")
for ax in axs:
    ax.set_xlim(0, 50); ax.set_xlabel("Tiempo de producción (ns)")
guardar(fig, "fig50_dssp")
print("listo:", sorted(os.listdir(FIG)))
