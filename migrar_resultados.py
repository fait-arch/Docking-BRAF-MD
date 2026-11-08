#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Agrupa las corridas sueltas de resultados/ en carpetas por sesión.

    python3 migrar_resultados.py            # enseña el plan, no mueve nada
    python3 migrar_resultados.py --aplicar  # lo ejecuta

Pasa de esto:

    resultados/10nsDabrafenib/
    resultados/10nsQuinine/
    resultados/0.1nsQuinineProtonada_20261108-0439/

a esto:

    resultados/20261106-0107/10nsDabrafenib/
    resultados/20261106-0107/10nsQuinine/
    resultados/20261108-0439/0.1nsQuinineProtonada/

La fecha de cada corrida sale del `run_id` de su estado/manifest.json, que es
el momento real del envío. Si falta, se usa el sufijo del nombre y, en último
caso, la fecha de modificación de la carpeta.

No pisa nada: si el destino ya existe, salta esa carpeta y lo avisa.
"""
import argparse
import datetime
import json
import os
import re
import shutil
import sys

RAIZ = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(RAIZ, "resultados")
RE_SESION = re.compile(r"^\d{8}-\d{4}$")
RE_CORRIDA = re.compile(r"^([\d.]+ns.+?)(?:_(\d{8}-\d{4}))?$")
#: compuesto-pipeline-tiempo-AAAAMMDD-HHMMSS
RE_RUNID = re.compile(r"(\d{8})-(\d{6})$")


def sello(carpeta, sufijo):
    """Marca de sesión de una corrida: manifiesto > sufijo > mtime."""
    mani = os.path.join(carpeta, "estado", "manifest.json")
    try:
        rid = json.load(open(mani)).get("run_id", "")
        m = RE_RUNID.search(rid)
        if m:
            return f"{m.group(1)}-{m.group(2)[:4]}", "manifiesto"
    except (OSError, ValueError):
        pass
    if sufijo:
        return sufijo, "nombre"
    t = datetime.datetime.fromtimestamp(os.path.getmtime(carpeta))
    return t.strftime("%Y%m%d-%H%M"), "fecha del directorio"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aplicar", action="store_true",
                    help="mueve de verdad; sin esto solo enseña el plan")
    a = ap.parse_args()

    if not os.path.isdir(BASE):
        sys.exit(f"No existe {BASE}")

    plan, saltadas = [], []
    for nombre in sorted(os.listdir(BASE)):
        origen = os.path.join(BASE, nombre)
        if not os.path.isdir(origen) or RE_SESION.match(nombre):
            continue                      # ya es una carpeta de sesión
        m = RE_CORRIDA.match(nombre)
        if not m:
            continue
        limpio, sufijo = m.group(1), m.group(2)
        s, de = sello(origen, sufijo)
        destino = os.path.join(BASE, s, limpio)
        (saltadas if os.path.exists(destino) else plan).append(
            (origen, destino, s, limpio, de))

    if not plan and not saltadas:
        print("\n  Nada que mover: todo está ya agrupado por sesión.\n")
        return

    print()
    for _, _, s, limpio, de in plan:
        print(f"  {limpio:34s} ->  {s}/{limpio}      (fecha segun {de})")
    for _, _, s, limpio, _ in saltadas:
        print(f"  {limpio:34s} ->  SALTADA: {s}/{limpio} ya existe")

    if not a.aplicar:
        print(f"\n  {len(plan)} carpeta(s) por mover. Nada se ha tocado.")
        print("  Para hacerlo:  python3 migrar_resultados.py --aplicar\n")
        return

    for origen, destino, _, limpio, _ in plan:
        os.makedirs(os.path.dirname(destino), exist_ok=True)
        shutil.move(origen, destino)
        print(f"  movida  {limpio}")
    print(f"\n  {len(plan)} carpeta(s) movida(s).\n")


if __name__ == "__main__":
    main()
