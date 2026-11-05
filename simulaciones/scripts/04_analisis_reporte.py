#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
04_analisis_reporte.py — análisis completos de la prueba de MD y figuras para el reporte.

Genera en --fig (PDF vectoriales) y --out (CSV + resumen.json):
  termodinamica   T, densidad, energía potencial y volumen vs tiempo (todas las fases)
  rmsd            RMSD Cα proteína y ligando vs tiempo
  rmsf            RMSF Cα por residuo (producción) con regiones funcionales de BRAF
  rg              radio de giro de la proteína
  contactos       frecuencia de contacto ligando–residuo (producción)
  distancias      distancias ligando ↔ Lys483, Cys532 (bisagra), Asp594 (DFG)
  enterramiento   fracción de la superficie del ligando enterrada en la proteína
  vina            afinidades de las poses de docking
  ligando_2d      estructura 2D del ligando (SMILES de sistema/sistema_info.json)

Uso:  python 04_analisis_reporte.py --base .. --md md_prueba --fig ../../../../../Img/Pallidine_24543347 --out ../analisis
"""
import argparse, glob, json, os, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OFFSET = 447                     # resid tleap 1 → residuo 448 de BRAF
SUP, TXT, TXT2, GRID, BAND = "#ffffff", "#0b0b0b", "#52514e", "#e4e3df", "#f0efec"
C1, C2, C3, C4, C5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
REGIONES = [  # (inicio, fin, nombre corto) numeración BRAF
    (464, 471, "Lazo P"), (492, 507, "Hélice αC"), (529, 535, "Bisagra"),
    (574, 581, "Lazo cat."), (594, 596, "DFG"), (597, 623, "Lazo de activación"),
]
FASES = [("02_cal", "Calentamiento"), ("03_eq", "Equilibrado"), ("04_prod_*", "Producción")]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 9.5,
    "axes.edgecolor": "#c9c8c2", "axes.labelcolor": TXT2, "xtick.color": TXT2, "ytick.color": TXT2,
    "text.color": TXT, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "lines.linewidth": 1.4, "savefig.bbox": "tight", "savefig.pad_inches": 0.04,
    "pdf.fonttype": 42,
})


def guardar(fig, carpeta, nombre):
    fig.savefig(os.path.join(carpeta, nombre + ".pdf"), facecolor=SUP)
    fig.savefig(os.path.join(carpeta, nombre + ".png"), dpi=170, facecolor=SUP)
    plt.close(fig)


def bandas_fases(ax, limites, etiquetas=True):
    total = limites[-1][1] - limites[0][0]
    for i, (a, b, n) in enumerate(limites):
        if n != "Producción":
            ax.axvspan(a, b, color=BAND, lw=0, zorder=0)
        if i:
            ax.axvline(a, color="#c9c8c2", lw=0.8, zorder=1)
    if etiquetas:
        import matplotlib.transforms as mt
        tr = mt.blended_transform_factory(ax.transData, ax.transAxes)
        corto = {"Calentamiento": "Calent.", "Equilibrado": "Equil.", "Producción": "Producción"}
        for a, b, n in limites:
            if (b - a) / total < 0.07:   # fase demasiado corta para rotularla
                continue
            ax.text((a + b) / 2, 1.02, corto.get(n, n), ha="center", va="bottom", fontsize=7.5,
                    color=TXT2, transform=tr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="..")
    ap.add_argument("--md", default="md_prueba", help="carpeta de la dinámica dentro de --base")
    ap.add_argument("--fig", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.fig, exist_ok=True); os.makedirs(a.out, exist_ok=True)
    S = os.path.join(a.base, "sistema"); M = os.path.join(a.base, a.md)
    R = {}

    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
    from MDAnalysis.lib.formats.libdcd import DCDFile

    top_full = os.path.join(S, "sistema.prmtop"); top_dry = os.path.join(S, "complejo_seco.prmtop")
    ref_full = mda.Universe(top_full, os.path.join(M, "01_min.pdb"), format="PDB")
    n_dry = mda.Universe(top_dry).atoms.n_atoms

    # ------------------------------------------------------------- termodinámica
    filas, t0 = [], 0.0
    for pref, fase in FASES:
        logs = sorted(glob.glob(os.path.join(M, pref + ".log")))
        for lf in logs:
            d = pd.read_csv(lf, sep="\t")
            d.columns = ["paso", "t", "Epot", "Ecin", "T", "V", "rho", "vel"][:len(d.columns)]
            d["fase"] = fase
            filas.append(d)
    th = pd.concat(filas, ignore_index=True)
    # tiempo continuo: el reloj sigue de calentamiento a equilibrado; la producción empieza en 0
    T_CAL = float(th[th.fase == "Calentamiento"].t.max())
    T_EQ = float(th[th.fase == "Equilibrado"].t.max())
    th["t_ps"] = np.where(th.fase == "Producción", th.t + T_EQ, th.t)
    T_FIN = float(th.t_ps.max())
    ESC, UNI = (1000.0, "ns") if T_FIN > 5000 else (1.0, "ps")
    R["duraciones_ps"] = {"calentamiento": T_CAL, "equilibrado": T_EQ - T_CAL, "produccion": T_FIN - T_EQ}
    th.to_csv(os.path.join(a.out, "termodinamica.csv"), index=False)
    lim = [(0, T_CAL / ESC, "Calentamiento"), (T_CAL / ESC, T_EQ / ESC, "Equilibrado"),
           (T_EQ / ESC, T_FIN / ESC, "Producción")]
    prod = th[th.fase == "Producción"]
    R["T_prod"] = [round(prod["T"].mean(), 1), round(prod["T"].std(), 1)]
    R["rho_prod"] = [round(prod.rho.mean(), 4), round(prod.rho.std(), 4)]
    R["Epot_prod_kJ"] = [round(prod.Epot.mean()), round(prod.Epot.std())]
    R["V_prod_nm3"] = [round(prod.V.mean(), 1), round(prod.V.std(), 1)]
    R["rho_inicio_cal"] = round(th.rho.iloc[0], 4)
    R["vel_ns_dia"] = {f: round(th[th.fase == f].vel.replace(0, np.nan).mean(), 2) for f in th.fase.unique()}

    fig, axs = plt.subplots(4, 1, figsize=(6.6, 7.6), sharex=True, gridspec_kw={"hspace": 0.36})
    specs = [("T", "Temperatura (K)", C1, 310, "objetivo 310 K"),
             ("rho", "Densidad (g/mL)", C3, None, None),
             ("Epot", "E. potencial (10⁵ kJ/mol)", C2, None, None),
             ("V", "Volumen (nm³)", C5, None, None)]
    for ax, (col, lab, c, ref, reflab) in zip(axs, specs):
        y = th[col] / 1e5 if col == "Epot" else th[col]
        ax.plot(th.t_ps / ESC, y, color=c, lw=1.3)
        ax.set_ylabel(lab, fontsize=8)
        if ref:
            ax.axhline(ref, color=TXT2, lw=0.8, ls=(0, (3, 3)))
            ax.text(T_FIN / ESC * 0.99, ref - 12, reflab, ha="right", va="top", fontsize=7.5, color=TXT2)
        bandas_fases(ax, lim, etiquetas=(ax is axs[0]))
    axs[-1].set_xlabel(f"Tiempo ({UNI})"); axs[-1].set_xlim(0, T_FIN / ESC)
    fig.align_ylabels(axs)
    guardar(fig, a.fig, "fig_termodinamica")

    # ------------------------------------------------------------- trayectorias
    def universo(pref):
        dcds = [d for d in sorted(glob.glob(os.path.join(M, pref + ".dcd"))) if os.path.getsize(d) > 0]
        with DCDFile(dcds[0]) as fh:
            n = fh.header["natoms"]
        return mda.Universe(top_full if n == ref_full.atoms.n_atoms else top_dry, dcds), n

    serie = []
    t_off = {"Calentamiento": 0.0, "Equilibrado": T_CAL, "Producción": T_EQ}
    ref_dry_atoms = ref_full.atoms[:n_dry]
    for pref, fase in FASES:
        u, n = universo(pref)
        refsel = ref_full.atoms if n == ref_full.atoms.n_atoms else ref_dry_atoms
        rr = rms.RMSD(u, refsel, select="protein and name CA",
                      groupselections=["resname LIG and not element H"]).run().results.rmsd
        prot = u.select_atoms("protein"); lig = u.select_atoms("resname LIG and not element H")
        k483 = u.select_atoms(f"resid {483-OFFSET} and name NZ")
        c532 = u.select_atoms(f"resid {532-OFFSET} and name N")
        d594 = u.select_atoms(f"resid {594-OFFSET} and name CG")
        for i, ts in enumerate(u.trajectory):
            lc = lig.center_of_mass()
            dmin = lambda g: float(np.min(np.linalg.norm(lig.positions[:, None] - g.positions[None], axis=2)))
            serie.append({"fase": fase, "t_ps": t_off[fase] + (i + 1) * u.trajectory.dt,
                          "rmsd_CA": rr[i, 2], "rmsd_lig": rr[i, 3], "rg": prot.radius_of_gyration(),
                          "d_K483": dmin(k483), "d_C532": dmin(c532), "d_D594": dmin(d594)})
    se = pd.DataFrame(serie)
    se = pd.concat([pd.DataFrame([{"fase": "Minimizada", "t_ps": 0.0, "rmsd_CA": 0.0, "rmsd_lig": 0.0}]), se],
                   ignore_index=True)
    se.to_csv(os.path.join(a.out, "series_estructurales.csv"), index=False)
    se["t_plot"] = se.t_ps / ESC
    sp = se[se.fase == "Producción"]
    for c in ("rmsd_CA", "rmsd_lig", "rg", "d_K483", "d_C532", "d_D594"):
        R[c + "_prod"] = [round(sp[c].mean(), 2), round(sp[c].std(), 2), round(sp[c].min(), 2), round(sp[c].max(), 2)]
    # deriva del RMSD Cα en producción (Å/ns) por regresión lineal
    pend = np.polyfit(sp.t_ps / 1000, sp.rmsd_CA, 1)[0]
    R["rmsd_CA_pendiente_A_por_ns"] = round(float(pend), 2)
    R["rmsd_CA_final"] = round(float(sp.rmsd_CA.iloc[-1]), 2)
    eq = se[se.fase == "Equilibrado"]; R["rmsd_lig_eq"] = [round(eq.rmsd_lig.mean(), 2), round(eq.rmsd_lig.max(), 2)]

    # RMSD
    fig, axs = plt.subplots(2, 1, figsize=(6.6, 4.6), sharex=True, gridspec_kw={"hspace": 0.38})
    for ax, col, c, tit in [(axs[0], "rmsd_CA", C1, "Proteína: carbonos α"),
                            (axs[1], "rmsd_lig", C2, "Ligando: átomos pesados")]:
        ax.plot(se.t_plot, se[col], color=c)
        ax.set_ylim(0, max(se[col].max() * 1.2, 0.5)); ax.set_ylabel("RMSD (Å)")
        ax.set_title(tit, loc="left", fontweight="bold")
        ax.set_title(f"producción: {sp[col].mean():.2f} ± {sp[col].std():.2f} Å", loc="right", color=TXT2, fontsize=8)
        bandas_fases(ax, lim, etiquetas=False)
    for ax in axs:
        for a0, b0, n in lim:
            if (b0 - a0) / (T_FIN / ESC) < 0.07:
                continue
            n = {"Calentamiento": "Calent.", "Equilibrado": "Equil."}.get(n, n)
            ax.text((a0 + b0) / 2, ax.get_ylim()[1] * 0.93, n, ha="center", va="top", fontsize=7.5, color=TXT2)
    axs[1].set_xlabel(f"Tiempo ({UNI})"); axs[1].set_xlim(0, T_FIN / ESC)
    guardar(fig, a.fig, "fig_rmsd")

    # Rg
    fig, ax = plt.subplots(figsize=(6.6, 2.3))
    s2 = se.dropna(subset=["rg"])
    ax.plot(s2.t_plot, s2.rg, color=C3)
    ax.set_ylabel("Radio de giro (Å)"); ax.set_xlabel(f"Tiempo ({UNI})"); ax.set_xlim(0, T_FIN / ESC)
    lo, hi = s2.rg.min(), s2.rg.max(); ax.set_ylim(lo - 0.4, hi + 0.4)
    bandas_fases(ax, lim)
    guardar(fig, a.fig, "fig_rg")

    # Distancias
    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    for col, c, lab in [("d_K483", C1, "Lys483 (NZ)"), ("d_C532", C2, "Cys532 (N, bisagra)"),
                        ("d_D594", C3, "Asp594 (CG, DFG)")]:
        ax.plot(se.t_plot, se[col], color=c, label=lab, lw=1.0 if len(se) > 500 else 1.4)
    ax.set_ylabel("Distancia mínima al ligando (Å)"); ax.set_xlabel(f"Tiempo ({UNI})"); ax.set_xlim(0, T_FIN / ESC)
    ax.set_ylim(0, se[["d_K483", "d_C532", "d_D594"]].max().max() * 1.15)
    bandas_fases(ax, lim)
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.22), fontsize=8)
    guardar(fig, a.fig, "fig_distancias")

    # ------------------------------------------------------------- RMSF y contactos (producción)
    u, _ = universo("04_prod_*")
    align.AlignTraj(u, u, select="protein and name CA", in_memory=True, ref_frame=0).run()
    ca = u.select_atoms("protein and name CA")
    # superponer sobre la estructura promedio
    from MDAnalysis.analysis.align import AverageStructure
    avg = AverageStructure(u, u, select="protein and name CA", ref_frame=0).run().results.universe
    align.AlignTraj(u, avg, select="protein and name CA", in_memory=True).run()
    rf = rms.RMSF(ca).run().results.rmsf
    resn = ca.resids + OFFSET
    rmsf = pd.DataFrame({"residuo": resn, "aa": ca.resnames, "rmsf_A": rf})
    rmsf.to_csv(os.path.join(a.out, "rmsf.csv"), index=False)
    top = rmsf.sort_values("rmsf_A", ascending=False).head(8)
    R["rmsf_media"] = round(float(rf.mean()), 2)
    R["rmsf_top"] = [(int(r), aa, round(float(v), 2)) for r, aa, v in zip(top.residuo, top.aa, top.rmsf_A)]
    R["rmsf_regiones"] = {n: round(float(rmsf[(rmsf.residuo >= i) & (rmsf.residuo <= f)].rmsf_A.mean()), 2)
                          for i, f, n in REGIONES}

    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    ymax = rf.max() * 1.25
    import matplotlib.transforms as mt
    tr = mt.blended_transform_factory(ax.transData, ax.transAxes)
    filas = [0, 1, 0, 1, 2, 0]
    for k, (i, f, n) in enumerate(REGIONES):
        ax.axvspan(i - 0.5, f + 0.5, color=BAND, lw=0, zorder=0)
        ax.text((i + f) / 2, 1.02 + 0.08 * filas[k], n, ha="center", va="bottom", fontsize=7,
                color=TXT2, transform=tr)
    ax.plot(resn, rf, color=C1, lw=1.2)
    ax.set_xlim(resn.min(), resn.max()); ax.set_ylim(0, ymax)
    ax.set_xlabel("Residuo de BRAF"); ax.set_ylabel("RMSF Cα (Å)")
    guardar(fig, a.fig, "fig_rmsf")

    # contactos
    from MDAnalysis.lib.distances import distance_array
    lig = u.select_atoms("resname LIG and not element H")
    prot = u.select_atoms("protein and not element H")
    cont = {}
    for ts in u.trajectory:
        d = distance_array(lig.positions, prot.positions)
        near = prot[np.any(d < 4.0, axis=0)]
        for r in set(zip(near.resids, near.resnames)):
            cont[r] = cont.get(r, 0) + 1
    nf = len(u.trajectory)
    cf = pd.DataFrame([{"residuo": r + OFFSET, "aa": n, "frecuencia": c / nf} for (r, n), c in cont.items()])
    cf = cf.sort_values("frecuencia", ascending=False)
    cf.to_csv(os.path.join(a.out, "contactos.csv"), index=False)
    R["contactos_50"] = [(int(r), aa, round(float(f), 2)) for r, aa, f in
                         zip(cf.residuo, cf.aa, cf.frecuencia) if f >= 0.5]

    # puentes de hidrógeno ligando-proteína
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
    hb_rows = []
    for d_sel, a_sel in [("resname LIG", "protein"), ("protein", "resname LIG")]:
        h = HBA(u, donors_sel=None, hydrogens_sel=None, acceptors_sel=None,
                between=["resname LIG", "protein"], d_a_cutoff=3.5, d_h_a_angle_cutoff=130)
        h.hydrogens_sel = h.guess_hydrogens("resname LIG or protein")
        h.acceptors_sel = h.guess_acceptors("resname LIG or protein")
        h.run()
        break
    hb = h.results.hbonds
    occ = {}
    for row in hb:
        d, hh, ac = int(row[1]), int(row[2]), int(row[3])
        key = (u.atoms[d].resname + str(u.atoms[d].resid + (0 if u.atoms[d].resname == "LIG" else OFFSET)) + ":" + u.atoms[d].name,
               u.atoms[ac].resname + str(u.atoms[ac].resid + (0 if u.atoms[ac].resname == "LIG" else OFFSET)) + ":" + u.atoms[ac].name)
        occ.setdefault(key, set()).add(int(row[0]))
    hbt = sorted([(k[0], k[1], len(v) / nf) for k, v in occ.items()], key=lambda x: -x[2])
    pd.DataFrame(hbt, columns=["donador", "aceptor", "ocupacion"]).to_csv(os.path.join(a.out, "puentes_H.csv"), index=False)
    R["puentes_H"] = [(d, ac, round(o, 2)) for d, ac, o in hbt[:6]]
    R["frames_con_puenteH"] = round(len(set(int(r[0]) for r in hb)) / nf, 2)

    fig, ax = plt.subplots(figsize=(6.6, 3.0))
    cfp = cf[cf.frecuencia >= 0.10].copy()
    nombre_aa = lambda x: {"HID": "His", "HIE": "His", "HIP": "His", "GLH": "Glu", "ASH": "Asp",
                           "CYX": "Cys", "LYN": "Lys"}.get(x, x.capitalize())
    etiquetas = [f"{nombre_aa(aa)}{r}" for aa, r in zip(cfp.aa, cfp.residuo)]
    colores = [C2 if f >= 0.5 else "#f2b199" for f in cfp.frecuencia]
    ax.bar(range(len(cfp)), cfp.frecuencia * 100, color=colores, width=0.72)
    ax.set_xticks(range(len(cfp))); ax.set_xticklabels(etiquetas, rotation=60, ha="right", fontsize=7.5)
    ax.set_ylabel("% de fotogramas en contacto"); ax.set_ylim(0, 105)
    ax.axhline(50, color=TXT2, lw=0.8, ls=(0, (3, 3)))
    ax.text(len(cfp) - 0.5, 51, "50 %", ha="right", va="bottom", fontsize=7.5, color=TXT2)
    ax.grid(axis="x", visible=False)
    guardar(fig, a.fig, "fig_contactos")

    # enterramiento del ligando (SASA)
    import mdtraj as md
    prod_dcds = sorted(glob.glob(os.path.join(M, "04_prod_*.dcd")))
    tj = md.load(prod_dcds, top=md.load_prmtop(top_dry))
    sel_c = tj.topology.select("protein or resname LIG"); tj = tj.atom_slice(sel_c)
    lig_idx = tj.topology.select("resname LIG")
    sasa_c = md.shrake_rupley(tj, mode="atom")[:, lig_idx].sum(1)
    sasa_l = md.shrake_rupley(tj.atom_slice(lig_idx), mode="atom").sum(1)
    ent = 1 - sasa_c / sasa_l
    R["ligando_enterrado_pct"] = [round(float(ent.mean() * 100), 1), round(float(ent.std() * 100), 1)]
    R["ligando_sasa_libre_A2"] = round(float(sasa_l.mean() * 100), 0)
    # estructura secundaria inicio vs final de la producción
    ss = md.compute_dssp(tj.atom_slice(tj.topology.select("protein"))[[0, -1]], simplified=True)
    R["dssp_inicio"] = {k: round(float((ss[0] == k).mean() * 100), 1) for k in ("H", "E", "C")}
    R["dssp_final"] = {k: round(float((ss[1] == k).mean() * 100), 1) for k in ("H", "E", "C")}

    # ------------------------------------------------------------- docking
    vals = []
    for l in open(os.path.join(S, "vina_poses.pdbqt")):
        if l.startswith("REMARK VINA RESULT"):
            p = l.split(); vals.append((float(p[3]), float(p[4]), float(p[5])))
    vd = pd.DataFrame(vals, columns=["afinidad", "rmsd_lb", "rmsd_ub"]); vd.index += 1
    vd.to_csv(os.path.join(a.out, "vina_poses.csv"), index_label="pose")
    R["vina"] = {"mejor": vd.afinidad.iloc[0], "n": len(vd), "top5": vd.afinidad.head(5).tolist(),
                 "peor": vd.afinidad.iloc[-1]}
    fig, ax = plt.subplots(figsize=(6.6, 2.4))
    ax.plot(vd.index, vd.afinidad, color=C1, marker="o", ms=2.8, lw=1.0)
    ax.scatter([1], [vd.afinidad.iloc[0]], s=40, color=C2, zorder=4)
    ax.annotate(f"pose usada: {vd.afinidad.iloc[0]:.2f} kcal/mol", (1, vd.afinidad.iloc[0]), xytext=(8, -4),
                textcoords="offset points", fontsize=8, color=TXT)
    ax.set_xlabel("Pose (ordenada por afinidad)"); ax.set_ylabel("Afinidad Vina (kcal/mol)")
    ax.set_xlim(0, len(vd) + 1)
    guardar(fig, a.fig, "fig_vina")

    # ------------------------------------------------------------- ligando 2D
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem, Descriptors, rdMolDescriptors, Crippen
    smi = json.load(open(os.path.join(S, "sistema_info.json")))["ligando_smiles"]
    m = Chem.MolFromSmiles(smi)
    from rdkit.Chem import rdDepictor
    rdDepictor.SetPreferCoordGen(True)   # dibuja mejor los anillos puenteados (morfinano)
    rdDepictor.Compute2DCoords(m)
    d2 = Draw.rdMolDraw2D.MolDraw2DCairo(900, 600)
    o = d2.drawOptions(); o.bondLineWidth = 3; o.minFontSize = 26; o.padding = 0.08
    d2.DrawMolecule(m); d2.FinishDrawing()
    open(os.path.join(a.fig, "fig_ligando_2d.png"), "wb").write(d2.GetDrawingText())
    R["ligando"] = {"formula": rdMolDescriptors.CalcMolFormula(m), "MW": round(Descriptors.MolWt(m), 2),
                    "logP": round(Crippen.MolLogP(m), 2), "HBD": rdMolDescriptors.CalcNumHBD(m),
                    "HBA": rdMolDescriptors.CalcNumHBA(m), "rotables": rdMolDescriptors.CalcNumRotatableBonds(m),
                    "TPSA": round(rdMolDescriptors.CalcTPSA(m), 1)}

    json.dump(R, open(os.path.join(a.out, "resumen.json"), "w"), indent=2, ensure_ascii=False, default=float)
    print(json.dumps(R, indent=2, ensure_ascii=False, default=float))


if __name__ == "__main__":
    main()
