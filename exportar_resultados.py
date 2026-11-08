#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reúne los resultados de las corridas en algo exportable.

    python3 exportar_resultados.py [tiempo]        # por defecto: 10

Produce dos cosas en la raíz del proyecto:

    exportacion/resumen_<tiempo>ns.csv   una fila por compuesto, todas las
                                         magnitudes comparables en columnas
    exportacion/<tiempo>ns.tar.gz        figuras, CSV y MM/GBSA de cada
                                         compuesto, SIN las trayectorias

Las trayectorias (md/) se excluyen a propósito: son el 99 % del peso y no
sirven para un informe. Si las necesitas, cópialas aparte con rsync.

No usa nada fuera de la librería estándar: no hace falta activar el entorno.
"""
import csv
import glob
import json
import os
import re
import sys
import tarfile

RAIZ = os.path.dirname(os.path.abspath(__file__))
TIEMPO = sys.argv[1] if len(sys.argv) > 1 else "10"
SALIDA = os.path.join(RAIZ, "exportacion")

#: Columnas del CSV.
CAMPOS_LIGANDO = ("formula", "MW", "logP", "HBD", "HBA", "TPSA", "rotables")
CAMPOS_RESUMEN = ("rmsd_CA_final", "rmsd_CA_pendiente_A_por_ns", "rmsf_media",
                  "frames_con_puenteH", "ligando_sasa_libre_A2")
#: Términos de MM/GBSA, en el bloque de diferencias (el último de cada archivo).
TERMINOS = ("VDWAALS", "EEL", "EGB", "ESURF")


def ultimo_valor(texto, termino):
    """El valor del término en el bloque 'Differences', que es el último."""
    hits = re.findall(rf"^\s*{termino}\s+(-?[\d.]+)", texto, re.M)
    return float(hits[-1]) if hits else None


def leer(carpeta):
    """Las magnitudes de una corrida, o None si está incompleta."""
    try:
        r = json.load(open(os.path.join(carpeta, "analisis", "resumen.json")))
    except (OSError, ValueError):
        return None
    fila = {"compuesto": os.path.basename(carpeta)[len(TIEMPO) + 2:]}
    lig = r.get("ligando", {})
    fila.update({c: lig.get(c) for c in CAMPOS_LIGANDO})
    fila["vina_mejor"] = r.get("vina", {}).get("mejor")
    fila["vina_peor"] = r.get("vina", {}).get("peor")
    fila.update({c: r.get(c) for c in CAMPOS_RESUMEN})

    dat = os.path.join(carpeta, "analisis", "mmgbsa", "FINAL_MMGBSA.dat")
    if os.path.exists(dat):
        texto = open(dat).read()
        for t in TERMINOS:
            fila[f"mmgbsa_{t.lower()}"] = ultimo_valor(texto, t)
        m = re.search(r"DELTA TOTAL\s+(-?[\d.]+)\s+([\d.]+)", texto)
        if m:
            fila["mmgbsa_total"] = float(m.group(1))
            fila["mmgbsa_sd"] = float(m.group(2))
        # La interacción directa, sin el término de solvatación, que es la
        # parte menos fiable del método.
        v, e = fila.get("mmgbsa_vdwaals"), fila.get("mmgbsa_eel")
        if v is not None and e is not None:
            fila["interaccion_vdw_mas_eel"] = round(v + e, 4)
    return fila


def main():
    carpetas = sorted(glob.glob(os.path.join(RAIZ, "resultados", f"{TIEMPO}ns*")))
    if not carpetas:
        sys.exit(f"No hay carpetas resultados/{TIEMPO}ns* en {RAIZ}")

    filas, incompletas = [], []
    for c in carpetas:
        f = leer(c)
        (filas if f else incompletas).append(f or os.path.basename(c))

    os.makedirs(SALIDA, exist_ok=True)

    columnas = ["compuesto", *CAMPOS_LIGANDO, "vina_mejor", "vina_peor",
                "mmgbsa_total", "mmgbsa_sd", "mmgbsa_vdwaals", "mmgbsa_eel",
                "mmgbsa_egb", "mmgbsa_esurf", "interaccion_vdw_mas_eel",
                *CAMPOS_RESUMEN]
    filas.sort(key=lambda f: f.get("mmgbsa_total") or 0)
    csv_ruta = os.path.join(SALIDA, f"resumen_{TIEMPO}ns.csv")
    with open(csv_ruta, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columnas, extrasaction="ignore")
        w.writeheader()
        w.writerows(filas)

    tar_ruta = os.path.join(SALIDA, f"{TIEMPO}ns.tar.gz")
    with tarfile.open(tar_ruta, "w:gz") as tar:
        for c in carpetas:
            for sub in ("resultados", "analisis", "estado"):
                d = os.path.join(c, sub)
                if os.path.isdir(d):
                    tar.add(d, arcname=os.path.join(os.path.basename(c), sub))
        tar.add(csv_ruta, arcname=os.path.basename(csv_ruta))

    print(f"\n  {'compuesto':20s} {'Vina':>7s} {'MM/GBSA':>16s} "
          f"{'vdW+EEL':>8s} {'EGB':>7s} {'RMSD':>5s} {'H-b':>5s}")
    print("  " + "-" * 76)
    for f in filas:
        print(f"  {f['compuesto']:20s} {f.get('vina_mejor') or 0:7.2f} "
              f"{f.get('mmgbsa_total') or 0:8.2f} ± {f.get('mmgbsa_sd') or 0:5.2f} "
              f"{f.get('interaccion_vdw_mas_eel') or 0:8.2f} "
              f"{f.get('mmgbsa_egb') or 0:7.2f} "
              f"{f.get('rmsd_CA_final') or 0:5.2f} "
              f"{f.get('frames_con_puenteH') or 0:5.2f}")
    if incompletas:
        print(f"\n  incompletas (sin resumen.json): {', '.join(incompletas)}")
    print(f"\n  {csv_ruta}")
    print(f"  {tar_ruta}  ({os.path.getsize(tar_ruta) / 1e6:.1f} MB)\n")


if __name__ == "__main__":
    main()
