#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
benchmark.py — mide la velocidad (ns/día) del sistema con los parámetros del artículo,
para estimar cuánto tardará la producción completa.

  python benchmark.py --sistema ../sistema --estado ../md_prueba/01_min.xml --pasos 1000
"""
import argparse, time
import openmm as mm
from openmm import app, unit as u

ap = argparse.ArgumentParser()
ap.add_argument("--sistema", default="../sistema")
ap.add_argument("--estado", default="../md_prueba/01_min.xml")
ap.add_argument("--pasos", type=int, default=1000)
ap.add_argument("--plataforma", default="auto")
a = ap.parse_args()

prm = app.AmberPrmtopFile(f"{a.sistema}/sistema.prmtop")
s = prm.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.1 * u.nanometer, constraints=app.HBonds)
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
sim = app.Simulation(prm.topology, s, integ, plat)
sim.loadState(a.estado)
sim.context.setVelocitiesToTemperature(310 * u.kelvin)
sim.step(100)
t = time.time(); sim.step(a.pasos); dt = time.time() - t
nsdia = a.pasos * 2e-6 / (dt / 86400)
print(f"Plataforma {plat.getName()}: {a.pasos} pasos en {dt:.1f} s -> {nsdia:.2f} ns/día")
for ns in (0.1, 1, 300, 1000):
    print(f"  {ns:>6} ns  ≈ {ns / nsdia * 24:,.1f} h")
