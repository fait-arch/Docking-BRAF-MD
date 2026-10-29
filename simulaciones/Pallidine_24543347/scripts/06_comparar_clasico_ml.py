#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
06_comparar_clasico_ml.py — compara la dinámica clásica (ligando con GAFF) con la ML/MM (ligando con
ANI-2x) del mismo sistema, en la misma ventana de tiempo de producción (por ejemplo, los primeros 20 ns).

Calcula lo mismo para las dos trayectorias:
  - RMSD de los Cα y del ligando (respecto a su estructura minimizada) y radio de giro;
  - puentes de hidrógeno ligando–proteína (D–A ≤ 3.5 Å, D–H···A ≥ 130°): ocupación y número medio;
  - contactos por residuo (átomos pesados a menos de 4 Å);
  - torsiones del ligando alrededor de enlaces rotables (lo que más cambia al pasar de GAFF a ANI-2x);
  - MM/GBSA de los fotogramas de cada corrida (mismo modelo GB en las dos).

Uso (entorno md-braf, desde scripts/):
  python 06_comparar_clasico_ml.py --clasico ../md_50ns --ml ../md_ml_ani2x \
         --gb-clasico ../analisis_50ns/mmgbsa --gb-ml ../analisis_ml_20ns/mmgbsa --t-max 20 \
         --out ../comparacion_ml_20ns

Con una sola réplica por método, una diferencia pequeña puede deberse al azar de la trayectoria; el
resumen lo indica.
"""
import argparse, glob, io, json, os, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.lib.distances import distance_array, calc_dihedrals
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

OFFSET = 447
SUP, TXT, TXT2, GRID = "#ffffff", "#0b0b0b", "#52514e", "#e4e3df"
COL = {"Clásico (GAFF)": "#2a78d6", "ML/MM (ANI-2x)": "#eb6834"}
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.edgecolor": "#c9c8c2", "axes.labelcolor": TXT2,
    "xtick.color": TXT2, "ytick.color": TXT2, "text.color": TXT, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.axisbelow": True, "savefig.bbox": "tight", "savefig.pad_inches": 0.04, "pdf.fonttype": 42,
    "legend.frameon": False,
})
NOMBRE_AA = {"HID": "His", "HIE": "His", "HIP": "His", "GLH": "Glu", "ASH": "Asp", "CYX": "Cys", "LYN": "Lys"}


def guardar(fig, out, nombre):
    fig.savefig(os.path.join(out, nombre + ".pdf"), facecolor=SUP)
    fig.savefig(os.path.join(out, nombre + ".png"), dpi=170, facecolor=SUP)
    plt.close(fig)


def etiqueta(at):
    if at.resname == "LIG":
        return f"LIG:{at.name}"
    return f"{NOMBRE_AA.get(at.resname, at.resname.capitalize())}{at.resid + OFFSET}:{at.name}"


def torsiones_ligando(u):
    """Torsiones A-B-C-D (átomos pesados) alrededor de enlaces B-C rotables (que no están en un anillo)."""
    lig = u.select_atoms("resname LIG and not element H")
    ids = set(lig.indices)
    vec = {i: set() for i in ids}
    for b in u.bonds:
        i, j = b.atoms.indices
        if i in ids and j in ids:
            vec[i].add(j); vec[j].add(i)

    def conectados_sin(i, j):  # ¿sigue unido el grafo si se quita el enlace i-j?
        vist, pila = {i}, [i]
        while pila:
            k = pila.pop()
            for m in vec[k]:
                if (k, m) in ((i, j), (j, i)) or m in vist:
                    continue
                vist.add(m); pila.append(m)
        return j in vist

    tors = []
    for b in sorted({tuple(sorted((i, j))) for i in ids for j in vec[i]}):
        B, C = b
        if conectados_sin(B, C):      # enlace en un anillo
            continue
        nb = sorted(vec[B] - {C}); nc = sorted(vec[C] - {B})
        if not nb or not nc:          # extremo (p. ej. un metilo): la torsión pesada no existe
            continue
        A, D = nb[0], nc[0]
        tors.append((A, B, C, D))
    return tors


def gbsa(carpeta, dt_ps, t_max_ns):
    if not carpeta or not os.path.exists(os.path.join(carpeta, "energias_por_frame.csv")):
        return None
    txt = open(os.path.join(carpeta, "energias_por_frame.csv")).read().splitlines()
    i0 = txt.index("DELTA Energy Terms") + 1
    filas = []
    for l in txt[i0:]:
        if not l.strip():
            break
        filas.append(l)
    E = pd.read_csv(io.StringIO("\n".join(filas)))
    intervalo, inicio = 1, 1
    for l in open(os.path.join(carpeta, "mmgbsa.in")):
        for p in l.replace(",", " ").split():
            if p.startswith("interval="):
                intervalo = int(p.split("=")[1])
            if p.startswith("startframe="):
                inicio = int(p.split("=")[1])
    E["t_ns"] = (inicio + intervalo * np.arange(len(E))) * dt_ps / 1000.0
    E = E[E.t_ns <= t_max_ns + 1e-6]
    return E.rename(columns={"VDWAALS": "vdw", "EEL": "eel", "EGB": "egb", "ESURF": "esurf",
                             "DELTA TOTAL": "dG"})[["t_ns", "dG", "vdw", "eel", "egb", "esurf"]]


def analizar(nombre, carpeta_md, top, t_max_ns, gb_dir):
    dcds = sorted(glob.glob(os.path.join(carpeta_md, "04_prod_*.dcd")))
    if not dcds:
        raise SystemExit(f"No hay trayectorias de producción en {carpeta_md}")
    u = mda.Universe(top, dcds)
    guardar_ps = u.trajectory.dt          # ps entre fotogramas (cabecera del DCD)
    nf = min(len(u.trajectory), int(round(t_max_ns * 1000 / guardar_ps)))
    if nf < 2:
        raise SystemExit(f"{carpeta_md}: con fotogramas cada {guardar_ps:g} ps, {t_max_ns:g} ns son menos de 2 fotogramas.")
    t_real = nf * guardar_ps / 1000
    ref_full = mda.Universe(os.path.join(carpeta_md, "01_min.pdb"))
    ref = ref_full.atoms[:u.atoms.n_atoms]
    rr = rms.RMSD(u, ref, select="protein and name CA",
                  groupselections=["resname LIG and not element H"]).run(stop=nf).results.rmsd
    prot = u.select_atoms("protein"); protH = u.select_atoms("protein and not element H")
    ligH = u.select_atoms("resname LIG and not element H")
    tors = torsiones_ligando(u)
    nombres_t = ["-".join(u.atoms[i].name for i in t) for t in tors]
    rg, cont, dih = [], {}, []
    for ts in u.trajectory[:nf]:
        rg.append(prot.radius_of_gyration())
        d = distance_array(ligH.positions, protH.positions)
        for r in set(protH[np.any(d < 4.0, axis=0)].residues):
            k = f"{NOMBRE_AA.get(r.resname, r.resname.capitalize())}{r.resid + OFFSET}"
            cont[k] = cont.get(k, 0) + 1
        if tors:
            p = u.atoms.positions
            dih.append(np.degrees(calc_dihedrals(*[p[[t[k] for t in tors]] for k in range(4)])))
    h = HBA(u, between=["resname LIG", "protein"], d_a_cutoff=3.5, d_h_a_angle_cutoff=130)
    h.hydrogens_sel = h.guess_hydrogens("resname LIG or protein")
    h.acceptors_sel = h.guess_acceptors("resname LIG or protein")
    h.run(stop=nf)
    hb = h.results.hbonds
    occ = {}
    for row in hb:
        occ.setdefault((etiqueta(u.atoms[int(row[1])]), etiqueta(u.atoms[int(row[3])])), set()).add(int(row[0]))
    S = pd.DataFrame({"t_ns": (np.arange(nf) + 1) * guardar_ps / 1000, "rmsd_CA": rr[:, 2], "rmsd_lig": rr[:, 3],
                      "rg": rg})
    nhb = np.zeros(nf)
    for row in hb:
        nhb[int(row[0])] += 1
    S["n_puentesH"] = nhb
    return {"nombre": nombre, "nf": nf, "t_ns": t_real, "serie": S,
            "contactos": {k: v / nf for k, v in cont.items()},
            "puentes": {k: len(v) / nf for k, v in occ.items()},
            "torsiones": (nombres_t, np.array(dih) if dih else np.zeros((0, 0))),
            "gb": gbsa(gb_dir, guardar_ps, t_max_ns), "dt_ps": guardar_ps}


def media_circular(x):
    r = np.radians(x)
    m = np.degrees(np.arctan2(np.sin(r).mean(), np.cos(r).mean()))
    R = np.hypot(np.sin(r).mean(), np.cos(r).mean())
    return m, np.degrees(np.sqrt(max(-2 * np.log(max(R, 1e-12)), 0)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--sistema", default="../sistema")
    ap.add_argument("--clasico", default="../md_50ns")
    ap.add_argument("--ml", default="../md_ml_ani2x")
    ap.add_argument("--gb-clasico", default="../analisis_50ns/mmgbsa")
    ap.add_argument("--gb-ml", default=None)
    ap.add_argument("--t-max", type=float, required=True, help="ns de producción a comparar")
    ap.add_argument("--out", required=True)
    ap.add_argument("--titulo", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    top = os.path.join(a.sistema, "complejo_seco.prmtop")
    R = [analizar("Clásico (GAFF)", a.clasico, top, a.t_max, a.gb_clasico),
         analizar("ML/MM (ANI-2x)", a.ml, top, a.t_max, a.gb_ml)]
    t_comun = min(r["t_ns"] for r in R)
    if t_comun < a.t_max - 1e-6:
        print(f"AVISO: una de las corridas solo tiene {t_comun:g} ns; se compara hasta ahí.")
        R = [analizar(r["nombre"], d, top, t_comun, g)
             for r, d, g in zip(R, (a.clasico, a.ml), (a.gb_clasico, a.gb_ml))]

    # ---------------------------------------------------------------- tabla resumen
    filas = []
    def fila(metrica, f, unidad=""):
        v = [f(r) for r in R]
        filas.append({"métrica": metrica, "unidad": unidad, R[0]["nombre"]: v[0], R[1]["nombre"]: v[1]})
    ms = lambda c: (lambda r: f"{r['serie'][c].mean():.2f} ± {r['serie'][c].std():.2f}")
    fila("RMSD Cα", ms("rmsd_CA"), "Å")
    fila("RMSD ligando", ms("rmsd_lig"), "Å")
    fila("Radio de giro", ms("rg"), "Å")
    fila("Puentes de H por fotograma", lambda r: f"{r['serie'].n_puentesH.mean():.2f}")
    fila("Fotogramas con algún puente de H", lambda r: f"{(r['serie'].n_puentesH > 0).mean() * 100:.0f} %")
    pares = sorted(set(R[0]["puentes"]) | set(R[1]["puentes"]),
                   key=lambda k: -max(R[0]["puentes"].get(k, 0), R[1]["puentes"].get(k, 0)))[:5]
    for k in pares:
        fila(f"Puente {k[0]} → {k[1]}", lambda r, k=k: f"{r['puentes'].get(k, 0) * 100:.0f} %")
    for r in R:
        r["n_cont50"] = sum(v >= 0.5 for v in r["contactos"].values())
    fila("Residuos en contacto > 50 % del tiempo", lambda r: f"{r['n_cont50']}")
    if all(r["gb"] is not None and len(r["gb"]) for r in R):
        for c, n in [("dG", "MM/GBSA ΔG"), ("vdw", "  van der Waals"), ("eel", "  electrostática"),
                     ("egb", "  solvatación polar (GB)")]:
            fila(n, lambda r, c=c: f"{r['gb'][c].mean():.1f} ± {r['gb'][c].std():.1f}", "kcal/mol")
        fila("  fotogramas de MM/GBSA", lambda r: f"{len(r['gb'])}")
    T = pd.DataFrame(filas)
    T.to_csv(os.path.join(a.out, "tabla_comparacion.csv"), index=False)

    # torsiones
    nt, d0 = R[0]["torsiones"]; _, d1 = R[1]["torsiones"]
    trows = []
    bins = np.linspace(-180, 180, 37)
    if d0.size == 0 or d1.size == 0:
        nt = []
    for j, n in enumerate(nt):
        p0, _ = np.histogram(d0[:, j], bins=bins); p1, _ = np.histogram(d1[:, j], bins=bins)
        sol = float(np.minimum(p0 / p0.sum(), p1 / p1.sum()).sum())
        m0, s0 = media_circular(d0[:, j]); m1, s1 = media_circular(d1[:, j])
        trows.append({"torsion": n, "media_clasico": round(m0, 1), "dispersion_clasico": round(s0, 1),
                      "media_ml": round(m1, 1), "dispersion_ml": round(s1, 1), "solapamiento": round(sol, 2)})
    TT = pd.DataFrame(trows)
    TT.to_csv(os.path.join(a.out, "torsiones_ligando.csv"), index=False)
    pd.DataFrame([{"donador": k[0], "aceptor": k[1], **{r["nombre"]: r["puentes"].get(k, 0) for r in R}}
                  for k in sorted(set(R[0]["puentes"]) | set(R[1]["puentes"]))]).to_csv(
        os.path.join(a.out, "puentes_H.csv"), index=False)
    C = pd.DataFrame({r["nombre"]: pd.Series(r["contactos"]) for r in R}).fillna(0)
    C = C.loc[C.max(axis=1).sort_values(ascending=False).index]
    C.to_csv(os.path.join(a.out, "contactos.csv"))

    # ---------------------------------------------------------------- figuras
    fig, axs = plt.subplots(3, 1, figsize=(7.0, 5.6), sharex=True, gridspec_kw={"hspace": 0.4})
    for r in R:
        s = r["serie"]; c = COL[r["nombre"]]
        for ax, col in zip(axs, ("rmsd_CA", "rmsd_lig", "n_puentesH")):
            y = s[col].rolling(50, center=True, min_periods=1).mean()
            ax.plot(s.t_ns, s[col], color=c, lw=0.4, alpha=0.35)
            ax.plot(s.t_ns, y, color=c, lw=1.4, label=r["nombre"])
    axs[0].set_ylabel("RMSD Cα (Å)"); axs[1].set_ylabel("RMSD ligando (Å)"); axs[2].set_ylabel("Puentes de H")
    axs[2].set_xlabel("Tiempo de producción (ns)"); axs[0].legend(loc="lower right", ncol=2)
    axs[0].set_title((a.titulo + " · " if a.titulo else "") + f"clásico frente a ML/MM, primeros {t_comun:g} ns",
                     loc="left", fontweight="bold")
    guardar(fig, a.out, "fig_cmp_estabilidad")

    top_c = C[C.max(axis=1) >= 0.25].head(16)
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.4), gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.65})
    x = np.arange(len(top_c)); w = 0.4
    for k, r in enumerate(R):
        axs[0].bar(x + (k - 0.5) * w, top_c[r["nombre"]] * 100, w, color=COL[r["nombre"]], label=r["nombre"])
    axs[0].set_xticks(x); axs[0].set_xticklabels(top_c.index, rotation=60, ha="right", fontsize=7.5)
    axs[0].set_ylabel("% de fotogramas en contacto"); axs[0].legend(loc="upper right", fontsize=7.5)
    axs[0].grid(axis="x", visible=False)
    if pares:
        y = np.arange(len(pares))
        for k, r in enumerate(R):
            axs[1].barh(y + (k - 0.5) * w, [r["puentes"].get(p, 0) * 100 for p in pares], w, color=COL[r["nombre"]])
        axs[1].set_yticks(y); axs[1].set_yticklabels([f"{p[0]}→{p[1]}" for p in pares], fontsize=7)
        axs[1].invert_yaxis(); axs[1].set_xlabel("Ocupación (%)"); axs[1].set_xlim(0, 100)
    axs[1].set_title("Puentes de hidrógeno", loc="left", fontsize=8.5)
    guardar(fig, a.out, "fig_cmp_interacciones")

    if len(nt):
        nc = min(4, len(nt)); nr = int(np.ceil(len(nt) / nc))
        fig, axs = plt.subplots(nr, nc, figsize=(1.85 * nc, 1.6 * nr + 0.3), squeeze=False)
        for j, n in enumerate(nt):
            ax = axs[j // nc, j % nc]
            for dd, r in ((d0, R[0]), (d1, R[1])):
                ax.hist(dd[:, j], bins=bins, density=True, histtype="stepfilled", alpha=0.45, color=COL[r["nombre"]])
            ax.set_title(f"{n}\nsolap. {TT.solapamiento[j]:.2f}", fontsize=6.8); ax.set_xlim(-180, 180)
            ax.set_xticks([-180, 0, 180]); ax.set_yticks([])
        for j in range(len(nt), nr * nc):
            axs[j // nc, j % nc].axis("off")
        fig.suptitle("Torsiones del ligando (azul: GAFF; naranja: ANI-2x)", fontsize=8.5)
        fig.tight_layout()
        guardar(fig, a.out, "fig_cmp_torsiones")

    if all(r["gb"] is not None and len(r["gb"]) for r in R):
        fig, ax = plt.subplots(figsize=(7.0, 2.6))
        for r in R:
            g = r["gb"]
            ax.plot(g.t_ns, g.dG, color=COL[r["nombre"]], lw=1.0,
                    label=f"{r['nombre']}: {g.dG.mean():.1f} ± {g.dG.std():.1f} kcal/mol")
        ax.set_xlabel("Tiempo de producción (ns)"); ax.set_ylabel("ΔG MM/GBSA (kcal/mol)"); ax.legend(loc="best")
        guardar(fig, a.out, "fig_cmp_mmgbsa")

    resumen = {"t_comparado_ns": t_comun, "fotogramas": [r["nf"] for r in R],
               "tabla": filas, "torsiones": trows,
               "nota": "Una réplica por método: diferencias pequeñas pueden deberse al azar de la trayectoria."}
    json.dump(resumen, open(os.path.join(a.out, "resumen_comparacion.json"), "w"), ensure_ascii=False, indent=1)
    pd.set_option("display.width", 160)
    print(f"\n== Comparación clásico vs ML/MM ({t_comun:g} ns de producción, {R[0]['nf']} fotogramas) ==")
    print(T.drop(columns="unidad").to_string(index=False))
    if len(TT):
        print("\nTorsiones del ligando (solapamiento 1 = misma distribución; < 0.5 = muy distinta):")
        print(TT.to_string(index=False))
    print("\nNota: una réplica por método; para afirmar una diferencia hacen falta varias réplicas.")
    print(f"Resultados en {a.out}")


if __name__ == "__main__":
    main()
