#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
benchmark.py — mide la velocidad (ns/día) del sistema con los parámetros del artículo,
para estimar cuánto tardará la producción completa.

  python benchmark.py --sistema ../sistema --estado ../md_prueba/01_min.xml --pasos 1000
  python benchmark.py --ml ani2x ...     # versión ML/MM (ligando con ANI-2x; entorno md-braf-ml)
"""
import argparse, time
import openmm as mm
from openmm import app, unit as u

ap = argparse.ArgumentParser()
ap.add_argument("--sistema", default="../sistema")
ap.add_argument("--estado", default="../md_prueba/01_min.xml")
ap.add_argument("--pasos", type=int, default=1000)
ap.add_argument("--plataforma", default="auto")
ap.add_argument("--ml", default="no", choices=("no", "ani2x", "ani1ccx"))
ap.add_argument("--ml-modelo-indice", type=int, default=-1)
ap.add_argument("--ml-doble", action="store_true")
ap.add_argument("--ml-cuda-graphs", action="store_true")
a = ap.parse_args()

prm = app.AmberPrmtopFile(f"{a.sistema}/sistema.prmtop")
s = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.1 * u.nanometer, constraints=app.HBonds)
if a.ml != "no":
    import tempfile, potencial_ml
    inp = app.AmberInpcrdFile(f"{a.sistema}/sistema.inpcrd")
    s, _ = potencial_ml.convertir_a_mlmm(s, prm.topology, inp.positions, a.ml, a.ml_modelo_indice,
                                         a.ml_cuda_graphs, tempfile.mkdtemp(), doble=a.ml_doble)
s.addForce(mm.MonteCarloBarostat(1 * u.atmosphere, 310 * u.kelvin, 25))
integ = mm.LangevinMiddleIntegrator(310 * u.kelvin, 1 / u.picosecond, 2 * u.femtosecond)
if a.plataforma == "auto":
    for n in ("CUDA", "HIP", "OpenCL", "CPU"):
        try:
            plat = mm.Platform.getPlatformByName(n); break
        except Exception:
            pass
else:
    plat = mm.Platform.getPlatformByName(a.plataforma)
props = {"Precision": "mixed"} if plat.getName() in ("CUDA", "OpenCL", "HIP") else {}
sim = app.Simulation(prm.topology, s, integ, plat, props)
sim.loadState(a.estado)
sim.context.setVelocitiesToTemperature(310 * u.kelvin)
sim.step(100)
t = time.time(); sim.step(a.pasos); dt = time.time() - t
nsdia = a.pasos * 2e-6 / (dt / 86400)
print(f"Plataforma {plat.getName()}{' + ' + a.ml.upper() if a.ml != 'no' else ''}: "
      f"{a.pasos} pasos en {dt:.1f} s -> {nsdia:.2f} ns/día")
for ns in (0.1, 1, 20, 50, 300, 1000):
    print(f"  {ns:>6} ns  ≈ {ns / nsdia * 24:,.1f} h")
