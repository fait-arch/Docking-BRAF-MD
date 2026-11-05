#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
05_preparar_mmgbsa.py — prepara la entrada de MM/GBSA para cualquier tiempo de
simulación.

Antes, el intervalo de fotogramas y el número de frames se calculaban en el
shell a partir del tiempo nominal (`T_PROD_PS / GUARDAR_PS`), lo que rompía el
cálculo si la producción se cortaba antes de terminar o si el tiempo cambiaba.
Aquí se cuentan los fotogramas que **realmente** hay en los .dcd y se deriva el
intervalo para quedarse con ~N fotogramas repartidos por toda la trayectoria.

Produce en --out:
  com.prmtop, rec.prmtop, lig.prmtop   (copiados si existen, o generados)
  mmgbsa.in                            (startframe/endframe/interval reales)
  dt_ps.txt                            (ps entre fotogramas usados, para la figura)
  frames.json                           (detalle de lo que se encontró)

Uso:
  python 05_preparar_mmgbsa.py --sistema ../sistema --md ../runs/50ns/md \\
         --out ../runs/50ns/analisis/mmgbsa --frames 200 --guardar-cada 10
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

PRMTOPS = ("com.prmtop", "rec.prmtop", "lig.prmtop")


def contar_fotogramas(top, dcds):
    """Número de fotogramas de la trayectoria concatenada."""
    try:
        import MDAnalysis as mda
    except ImportError:
        return None
    try:
        u = mda.Universe(top, list(dcds))
        return int(len(u.trajectory))
    except Exception as exc:                                    # noqa: BLE001
        print(f"  aviso: no pude leer la trayectoria con MDAnalysis ({exc})")
        return None


def localizar_prmtops(base, out):
    """Busca los prmtop de MM/GBSA ya generados en otra ejecución del compuesto."""
    candidatos = []
    candidatos += sorted(glob.glob(os.path.join(base, "analisis*", "mmgbsa")))
    candidatos += sorted(glob.glob(os.path.join(base, "runs", "*", "analisis", "mmgbsa")))
    for d in candidatos:
        if os.path.abspath(d) == os.path.abspath(out):
            continue
        if all(os.path.exists(os.path.join(d, p)) for p in PRMTOPS):
            return d
    return None


def generar_prmtops(sistema, out, ligando_resname):
    """Genera com/rec/lig.prmtop con ante-MMPBSA.py a partir del complejo seco."""
    com = os.path.join(sistema, "complejo_seco.prmtop")
    if not os.path.exists(com):
        sys.exit(f"ERROR: no existe {com} y tampoco encontré prmtop de MM/GBSA que copiar.\n"
                 "       Ejecuta antes la etapa 'sistema' del pipeline de resultados.")
    if shutil.which("ante-MMPBSA.py") is None:
        sys.exit("ERROR: falta ante-MMPBSA.py (AmberTools) para generar los prmtop.")
    shutil.copy(com, os.path.join(out, "com.prmtop"))
    cmd = ["ante-MMPBSA.py", "-p", com,
           "-c", os.path.join(out, "com.prmtop"),
           "-r", os.path.join(out, "rec.prmtop"),
           "-l", os.path.join(out, "lig.prmtop"),
           "-s", "", "-n", f":{ligando_resname}"]
    print("  " + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit("ERROR en ante-MMPBSA.py:\n" + (proc.stderr or proc.stdout))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--sistema", required=True)
    ap.add_argument("--md", required=True, help="carpeta con 04_prod_*.dcd")
    ap.add_argument("--out", required=True, help="carpeta de salida (analisis/mmgbsa)")
    ap.add_argument("--frames", type=int, default=200,
                    help="fotogramas objetivo para MM/GBSA")
    ap.add_argument("--guardar-cada", type=float, default=10.0,
                    help="ps entre fotogramas de la trayectoria")
    ap.add_argument("--ligando-resname", default="LIG")
    ap.add_argument("--titulo", default="MM/GBSA")
    ap.add_argument("--igb", type=int, default=5)
    ap.add_argument("--saltcon", type=float, default=0.150)
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    base = os.path.abspath(os.path.join(a.md, os.pardir))
    if os.path.basename(base) == "runs" or os.path.isdir(os.path.join(base, "sistema")) is False:
        # --md puede ser ../runs/50ns/md: el directorio del compuesto está 3 niveles arriba
        base = os.path.abspath(os.path.join(a.sistema, os.pardir))

    dcds = sorted(glob.glob(os.path.join(a.md, "04_prod_*.dcd")))
    if not dcds:
        sys.exit(f"ERROR: no hay trayectorias 04_prod_*.dcd en {a.md}")
    print(f"[1] Trayectorias: {len(dcds)} archivo(s)")

    top_seco = os.path.join(a.sistema, "complejo_seco.prmtop")
    n_frames = contar_fotogramas(top_seco, dcds)
    if n_frames is None:
        # reserva: estimar por el tiempo nominal de los logs
        n_frames = 0
        for d in dcds:
            log = d[:-4] + ".log"
            if os.path.exists(log):
                with open(log) as fh:
                    n_frames += max(sum(1 for _ in fh) - 1, 0)
        if n_frames == 0:
            sys.exit("ERROR: no pude determinar el número de fotogramas.")
        print(f"    (estimado por los .log: {n_frames})")
    print(f"    fotogramas disponibles: {n_frames}")

    intervalo = max(1, n_frames // max(a.frames, 1))
    usados = len(range(0, n_frames, intervalo))
    dt_ps = intervalo * a.guardar_cada
    print(f"[2] interval={intervalo} -> {usados} fotogramas, uno cada {dt_ps:g} ps")

    print("[3] Topologías de MM/GBSA")
    origen = localizar_prmtops(base, a.out)
    if origen:
        for p in PRMTOPS:
            shutil.copy(os.path.join(origen, p), os.path.join(a.out, p))
        print(f"    copiadas de {origen}")
    elif all(os.path.exists(os.path.join(a.out, p)) for p in PRMTOPS):
        print("    ya estaban en la carpeta de salida")
    else:
        print("    generando con ante-MMPBSA.py")
        generar_prmtops(a.sistema, a.out, a.ligando_resname)

    with open(os.path.join(a.out, "mmgbsa.in"), "w") as fh:
        fh.write(f"""{a.titulo}
&general
  startframe=1, endframe={n_frames}, interval={intervalo}, verbose=2, keep_files=0,
/
&gb
  igb={a.igb}, saltcon={a.saltcon},
/
&decomp
  idecomp=1, dec_verbose=0,
/
""")
    with open(os.path.join(a.out, "dt_ps.txt"), "w") as fh:
        fh.write(f"{dt_ps:g}\n")
    with open(os.path.join(a.out, "frames.json"), "w") as fh:
        json.dump({"trayectorias": [os.path.basename(d) for d in dcds],
                   "fotogramas_totales": n_frames,
                   "intervalo": intervalo,
                   "fotogramas_usados": usados,
                   "guardar_cada_ps": a.guardar_cada,
                   "dt_ps": dt_ps}, fh, indent=2)
    print(f"[4] Listo: {a.out}/mmgbsa.in")


if __name__ == "__main__":
    main()
