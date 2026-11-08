#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figura combinada de RMSD: los compuestos de un tiempo dado, superpuestos.

    python3 figura_rmsd.py [tiempo]          # por defecto: 10

Lee resultados/<tiempo>ns*/resultados/rmsd.csv y escribe en exportacion/figuras/
PNG (200 dpi) y PDF vectorial de:

    fig_rmsd_CA_combinado      el esqueleto de la proteína, las 5 superpuestas
    fig_rmsd_ligando_combinado el ligando dentro del bolsillo
    fig_rmsd_CA_paneles        lo mismo en paneles separados, eje compartido

Sobre los colores: la paleta está validada para daltonismo con todos los pares
en juego (peor par en visión normal ΔE 15.6, por encima del suelo de 15; CVD
6.9, dentro de la banda que exige codificación secundaria). Por eso cada serie
lleva ADEMÁS su propio patrón de trazo y una etiqueta directa al final: la
identidad nunca depende solo del color.
"""
import csv
import glob
import os
import sys

import matplotlib
import numpy as np
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RAIZ = os.path.dirname(os.path.abspath(__file__))
TIEMPO = sys.argv[1] if len(sys.argv) > 1 else "10"
SALIDA = os.path.join(RAIZ, "exportacion", "figuras")

#: Orden fijo: el color sigue al compuesto, nunca a su posición en el ranking.
#: El control positivo va primero.
ESTILO = {
    "Dabrafenib":      ("#2a78d6", (0, ()),              "Dabrafenib (control)"),
    "Macaridine":      ("#1baf7a", (0, (6, 2)),          "Macaridine"),
    "NotholaenicAcid": ("#eda100", (0, (1, 1.6)),        "Ácido notoloénico"),
    "Pallidine":       ("#008300", (0, (8, 2, 1.5, 2)),  "Pallidine"),
    "Quinine":         ("#e34948", (0, (3, 1.6)),        "Quinine"),
}
TINTA, TINTA2, REJILLA = "#1a1a19", "#5c5b55", "#e5e4df"


def series(carpeta, columna):
    """(tiempo en ns desde el inicio de producción, valor) de la fase prod."""
    t, v = [], []
    with open(os.path.join(carpeta, "resultados", "rmsd.csv")) as fh:
        for fila in csv.DictReader(fh):
            if fila["fase"] != "Producción":
                continue
            t.append((float(fila["tiempo_ps"]) - 1000.0) / 1000.0)
            v.append(float(fila[columna]))
    return t, v


def suavizar(v, ventana=51):
    """Media móvil centrada. Mil puntos crudos de cinco series son una maraña:
    la traza cruda va tenue al fondo y esta media al frente."""
    if len(v) < ventana:
        return np.asarray(v, float)
    k = np.ones(ventana) / ventana
    pad = ventana // 2
    ext = np.concatenate([np.full(pad, v[0]), np.asarray(v, float),
                          np.full(pad, v[-1])])
    return np.convolve(ext, k, mode="valid")


def estetica(ax, titulo, ylab, nota=None):
    ax.set_title(titulo, color=TINTA, fontsize=12, loc="left", pad=12)
    ax.set_xlabel("Tiempo de producción (ns)", color=TINTA2, fontsize=10)
    ax.set_ylabel(ylab, color=TINTA2, fontsize=10)
    ax.grid(axis="y", color=REJILLA, lw=0.8)
    ax.set_axisbelow(True)
    for lado in ("top", "right"):
        ax.spines[lado].set_visible(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_color(REJILLA)
    ax.tick_params(colors=TINTA2, labelsize=9, length=0)
    if nota:
        ax.text(0, -0.17, nota, transform=ax.transAxes, color=TINTA2,
                fontsize=8.5, va="top")


def etiquetas_sin_choque(ax, puntos, margen_rel=0.045):
    """Separa verticalmente las etiquetas directas que caerían encima."""
    y0, y1 = ax.get_ylim()
    sep = (y1 - y0) * margen_rel
    puntos.sort(key=lambda p: p[1])
    for i in range(1, len(puntos)):
        if puntos[i][1] - puntos[i - 1][1] < sep:
            puntos[i] = (puntos[i][0], puntos[i - 1][1] + sep, *puntos[i][2:])
    for x, y, color, texto in puntos:
        ax.plot(x + 0.12, y, "o", color=color, ms=4.5, clip_on=False, zorder=5)
        ax.text(x + 0.30, y, texto, color=TINTA, fontsize=8.5,
                va="center", clip_on=False)


def combinada(carpetas, columna, titulo, ylab, nota, nombre):
    fig, ax = plt.subplots(figsize=(9.4, 5.2))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    etiquetas, lo, hi = [], [], []
    for clave, (color, trazo, rotulo) in ESTILO.items():
        c = next((k for k in carpetas if k.endswith(clave)), None)
        if not c:
            continue
        t, v = series(c, columna)
        s = suavizar(v)
        ax.plot(t, v, color=color, lw=0.7, alpha=0.16, zorder=1)
        ax.plot(t, s, color=color, lw=2.2, linestyle=trazo, label=rotulo,
                solid_capstyle="round", zorder=3)
        etiquetas.append((t[-1], float(s[-1]), color, rotulo))
        lo.append(min(v)); hi.append(max(v))
    ax.set_xlim(0, float(TIEMPO))
    # Hueco arriba para que la leyenda no se apoye sobre los datos.
    y0, y1 = min(lo), max(hi)
    ax.set_ylim(y0 - 0.04 * (y1 - y0), y1 + 0.42 * (y1 - y0))
    estetica(ax, titulo, ylab, nota)
    leg = ax.legend(loc="upper left", frameon=False, fontsize=9,
                    ncol=3, handlelength=3.0, labelspacing=0.45,
                    columnspacing=1.6)
    for txt in leg.get_texts():
        txt.set_color(TINTA)
    etiquetas_sin_choque(ax, etiquetas)
    fig.subplots_adjust(left=0.08, right=0.80, top=0.90, bottom=0.17)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(SALIDA, f"{nombre}.{ext}"), dpi=200,
                    facecolor="white")
    plt.close(fig)


def paneles(carpetas, columna, titulo, ylab, nombre):
    presentes = [(k, v) for k, v in ESTILO.items()
                 if any(c.endswith(k) for c in carpetas)]
    fig, axes = plt.subplots(len(presentes), 1, figsize=(8.4, 1.5 * len(presentes)),
                             sharex=True, sharey=True)
    fig.patch.set_facecolor("white")
    for ax, (clave, (color, trazo, rotulo)) in zip(axes, presentes):
        c = next(k for k in carpetas if k.endswith(clave))
        t, v = series(c, columna)
        ax.plot(t, v, color=color, lw=0.6, alpha=0.20)
        ax.plot(t, suavizar(v), color=color, lw=1.8, linestyle=trazo)
        ax.set_facecolor("white")
        ax.grid(axis="y", color=REJILLA, lw=0.8)
        ax.set_axisbelow(True)
        for lado in ("top", "right"):
            ax.spines[lado].set_visible(False)
        for lado in ("left", "bottom"):
            ax.spines[lado].set_color(REJILLA)
        ax.tick_params(colors=TINTA2, labelsize=8.5, length=0)
        ax.text(0.012, 0.86, rotulo, transform=ax.transAxes, color=TINTA,
                fontsize=9, va="top")
    axes[0].set_title(titulo, color=TINTA, fontsize=12, loc="left", pad=12)
    axes[-1].set_xlabel("Tiempo de producción (ns)", color=TINTA2, fontsize=10)
    axes[len(axes) // 2].set_ylabel(ylab, color=TINTA2, fontsize=10)
    axes[0].set_xlim(0, float(TIEMPO))
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(SALIDA, f"{nombre}.{ext}"), dpi=200,
                    facecolor="white")
    plt.close(fig)


def main():
    carpetas = [c for c in sorted(glob.glob(
        os.path.join(RAIZ, "resultados", f"{TIEMPO}ns*"))) if os.path.isdir(c)]
    if not carpetas:
        sys.exit(f"No hay carpetas resultados/{TIEMPO}ns* en {RAIZ}")
    os.makedirs(SALIDA, exist_ok=True)

    nota = (f"BRAF V600E · {TIEMPO} ns de producción a 310 K y 1 atm · "
            "RMSD respecto a la estructura minimizada, tras alineamiento.")
    combinada(carpetas, "rmsd_CA_A",
              "Estabilidad del esqueleto proteico durante la dinámica",
              "RMSD de Cα (Å)", nota, "fig_rmsd_CA_combinado")
    combinada(carpetas, "rmsd_ligando_A",
              "Movilidad del ligando dentro del bolsillo de ATP",
              "RMSD del ligando (Å)", nota, "fig_rmsd_ligando_combinado")
    paneles(carpetas, "rmsd_CA_A",
            "RMSD de Cα, un panel por compuesto",
            "RMSD de Cα (Å)", "fig_rmsd_CA_paneles")

    for f in sorted(glob.glob(os.path.join(SALIDA, "*"))):
        print(f"  {os.path.getsize(f)/1024:7.0f} KB  {f}")


if __name__ == "__main__":
    main()
