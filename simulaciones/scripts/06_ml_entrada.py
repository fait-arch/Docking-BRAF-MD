#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
06_ml_entrada.py — puente entre la dinámica molecular y el pipeline de machine
learning.

Dos pasos, seleccionables con --paso:

  centro     calcula el centro de la caja de docking a partir de los residuos
             del bolsillo de ATP y lo escribe en centro_caja.txt ("x y z").
             Es lo que consume `dml-dock --center $CENTRO`.

  complejos  construye los PDB de complejo receptor+ligando que necesita
             `dml-features` y escribe complejos.txt (dos columnas:
             <pdb> <resname>). El origen depende de --fuente:

               --fuente md       N fotogramas repartidos por la producción del
                                 tiempo de simulación indicado (10, 20, 50 ns).
                                 Permite comparar el rescoring según la ventana
                                 simulada, sobre los mismos compuestos.
               --fuente docking  las poses acopladas en <out>/docked.

Uso:
  python 06_ml_entrada.py --paso centro   --sistema ../sistema --out ../runs/50ns/ml \\
         --bolsillo 463,471,481,483,505,514,529,531,532,583,594
  python 06_ml_entrada.py --paso complejos --fuente md --sistema ../sistema \\
         --md ../runs/50ns/md --out ../runs/50ns/ml --frames 20
"""
import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

RECEPTORES = ("receptor_amber.pdb", "receptor_H.pdb", "receptor_recortado.pdb")


# ---------------------------------------------------------------------------
#  centro del bolsillo
# ---------------------------------------------------------------------------
def leer_atomos(pdb):
    with open(pdb) as fh:
        for l in fh:
            if l.startswith(("ATOM", "HETATM")):
                yield l


def centro_bolsillo(pdb, residuos):
    xyz = []
    for l in leer_atomos(pdb):
        try:
            num = int(l[22:26])
        except ValueError:
            continue
        if num in residuos:
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    if not xyz:
        sys.exit(f"ERROR: ninguno de los residuos del bolsillo aparece en {pdb}")
    n = len(xyz)
    return [sum(p[i] for p in xyz) / n for i in range(3)]


def paso_centro(a):
    receptor = None
    for nombre in RECEPTORES:
        cand = os.path.join(a.sistema, nombre)
        if os.path.exists(cand):
            receptor = cand
            break
    if receptor is None:
        sys.exit(f"ERROR: no encontré el receptor preparado en {a.sistema} "
                 f"(busqué: {', '.join(RECEPTORES)})")
    residuos = {int(x) for x in a.bolsillo.split(",") if x.strip()}
    centro = centro_bolsillo(receptor, residuos)
    texto = " ".join(f"{c:.3f}" for c in centro)
    with open(os.path.join(a.out, "centro_caja.txt"), "w") as fh:
        fh.write(texto + "\n")
    with open(os.path.join(a.out, "caja.json"), "w") as fh:
        json.dump({"receptor": os.path.basename(receptor),
                   "residuos_bolsillo": sorted(residuos),
                   "centro": [round(c, 3) for c in centro]}, fh, indent=2)
    print(f"[centro] {receptor} -> centro = {texto}")


# ---------------------------------------------------------------------------
#  complejos a partir de la dinámica
# ---------------------------------------------------------------------------
def complejos_desde_md(a):
    try:
        import MDAnalysis as mda
    except ImportError:
        sys.exit("ERROR: MDAnalysis es necesario para --fuente md")

    dcds = sorted(glob.glob(os.path.join(a.md, "04_prod_*.dcd")))
    if not dcds:
        sys.exit(f"ERROR: no hay trayectorias 04_prod_*.dcd en {a.md}.\n"
                 "       Ejecuta antes el pipeline 'resultados' para este tiempo.")
    top = os.path.join(a.sistema, "complejo_seco.prmtop")
    if not os.path.exists(top):
        sys.exit(f"ERROR: falta {top}")

    u = mda.Universe(top, dcds)
    total = len(u.trajectory)
    n = max(1, min(a.frames, total))
    # fotogramas repartidos uniformemente por toda la producción
    indices = [round(i * (total - 1) / (n - 1)) if n > 1 else total - 1
               for i in range(n)]
    destino = os.path.join(a.out, "complejos")
    os.makedirs(destino, exist_ok=True)

    sel = u.select_atoms("all")
    salidas = []
    for k, idx in enumerate(indices, start=1):
        u.trajectory[int(idx)]
        pdb = os.path.join(destino, f"frame_{k:03d}.pdb")
        sel.write(pdb)
        salidas.append(pdb)
    print(f"[complejos] {len(salidas)} fotogramas de {total} "
          f"({os.path.basename(a.md)}) -> {destino}")
    return salidas, {"fuente": "md", "trayectorias": len(dcds),
                     "fotogramas_totales": total, "indices": indices}


# ---------------------------------------------------------------------------
#  complejos a partir del docking
# ---------------------------------------------------------------------------
def _concat_pdb(receptor, ligando, salida):
    with open(salida, "w") as out:
        for f in (receptor, ligando):
            with open(f) as fh:
                for l in fh:
                    if l.startswith(("ATOM", "HETATM")):
                        out.write(l)
            out.write("TER\n")
        out.write("END\n")


def complejos_desde_docking(a):
    if shutil.which("obabel") is None:
        sys.exit("ERROR: obabel es necesario para convertir las poses de docking")
    receptor = None
    for nombre in RECEPTORES:
        cand = os.path.join(a.sistema, nombre)
        if os.path.exists(cand):
            receptor = cand
            break
    if receptor is None:
        sys.exit(f"ERROR: no encontré el receptor preparado en {a.sistema}")

    poses = sorted(glob.glob(os.path.join(a.out, "docked", "*_out.pdbqt")))
    if not poses:
        poses = sorted(glob.glob(os.path.join(a.sistema, "vina_poses.pdbqt")))
    if not poses:
        sys.exit(f"ERROR: no hay poses acopladas en {a.out}/docked. "
                 "Ejecuta antes la etapa 'docking'.")

    destino = os.path.join(a.out, "complejos")
    os.makedirs(destino, exist_ok=True)
    salidas = []
    for pose_pdbqt in poses:
        etiqueta = os.path.basename(pose_pdbqt).replace("_out.pdbqt", "").replace(".pdbqt", "")
        tmp = os.path.join(destino, f"{etiqueta}_pose.pdb")
        # -m separa cada modo en su propio archivo: pose1.pdb, pose2.pdb, ...
        proc = subprocess.run(["obabel", pose_pdbqt, "-O", tmp, "-m"],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            print(f"  aviso: obabel falló con {pose_pdbqt}: {proc.stderr.strip()}")
            continue
        modos = sorted(glob.glob(os.path.join(destino, f"{etiqueta}_pose*.pdb")))
        for k, lig_pdb in enumerate(modos[:a.frames], start=1):
            comp = os.path.join(destino, f"{etiqueta}_m{k:02d}_complejo.pdb")
            _concat_pdb(receptor, lig_pdb, comp)
            salidas.append(comp)
    if not salidas:
        sys.exit("ERROR: no se pudo construir ningún complejo a partir del docking")
    print(f"[complejos] {len(salidas)} complejos de docking -> {destino}")
    return salidas, {"fuente": "docking", "poses": len(poses)}


def paso_complejos(a):
    if a.fuente == "md":
        salidas, meta = complejos_desde_md(a)
    else:
        salidas, meta = complejos_desde_docking(a)
    lista = os.path.join(a.out, "complejos.txt")
    with open(lista, "w") as fh:
        fh.write("# <pdb_complejo> <resname_ligando>\n")
        for p in salidas:
            fh.write(f"{p} {a.ligando_resname}\n")
    meta["ligando_resname"] = a.ligando_resname
    meta["complejos"] = len(salidas)
    with open(os.path.join(a.out, "complejos.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"[complejos] lista para dml-features: {lista}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--paso", default="centro", choices=["centro", "complejos"])
    ap.add_argument("--fuente", default="docking", choices=["docking", "md"])
    ap.add_argument("--sistema", required=True)
    ap.add_argument("--md", default=None, help="carpeta de la dinámica (--fuente md)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--frames", type=int, default=20,
                    help="fotogramas (md) o poses por ligando (docking)")
    ap.add_argument("--ligando-resname", default="LIG")
    ap.add_argument("--bolsillo",
                    default="463,471,481,483,505,514,529,531,532,583,594")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.paso == "centro":
        paso_centro(a)
    else:
        if a.fuente == "md" and not a.md:
            sys.exit("ERROR: --fuente md necesita --md")
        paso_complejos(a)


if __name__ == "__main__":
    main()
