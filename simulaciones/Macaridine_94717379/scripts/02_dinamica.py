#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
02_dinamica.py  —  protocolo minimización → calentamiento → equilibrado → producción
(RSC Adv. 2022, 12, 30181, sección 2.1.3), implementado en OpenMM sobre la
topología AMBER generada por 01_preparar_sistema.py.

Parámetros del artículo (valores por defecto de este script):
  - Minimización 10 000 ciclos con restricción armónica k = 10 kcal/mol/Å² sobre
    proteína + ligando (átomos pesados).  El artículo usa 4000 SD + 6000 CG en
    sander; aquí se usa el minimizador L-BFGS de OpenMM con el mismo límite de
    10 000 iteraciones y las mismas restricciones.
  - Calentamiento escalonado 0 → 310 K en 500 ps, sin restricciones (NVT).
  - Equilibrado 500 ps a 310 K y 1 atm (NPT).
  - Producción NPT 310 K, 1 atm (1 µs en el artículo).
  - Termostato de Langevin, frecuencia de colisión 1 ps⁻¹; paso de 2 fs;
    enlaces con H restringidos (equivalente a SHAKE); corte no enlazado 11 Å;
    electrostática de largo alcance PME; condiciones periódicas en todas direcciones.

Uso (protocolo completo del artículo — requiere GPU):
  python 02_dinamica.py --sistema ../sistema --out ../md

Prueba rápida (mismos parámetros, tiempos cortos):
  python 02_dinamica.py --sistema ../sistema --out ../md_prueba \
         --t-calentamiento 20 --t-equilibrado 20 --t-produccion 100 \
         --segmentos 4 --guardar-cada 1

Se puede ejecutar por etapas con --etapas min,cal,eq,prod (cada etapa parte del .xml
de la anterior). La producción se divide en segmentos 04_prod_001, 04_prod_002, ...
y es reanudable: si se interrumpe, al volver a lanzarla salta los segmentos ya hechos.

Versión con machine learning (ML/MM): --ml ani2x describe el ligando con el potencial de red
neuronal ANI-2x (ver potencial_ml.py) y deja todo lo demás igual. Necesita el entorno
md-braf-ml (environment-ml.yml) y un ligando neutro. Ejemplo:
  python 02_dinamica.py --ml ani2x --sistema ../sistema --out ../md_ml --t-produccion 20000 --segmentos 4
"""
import argparse, inspect, os, sys, time
import openmm as mm
from openmm import app, unit as u

KCAL_A2_TO_KJ_NM2 = 4.184 * 100.0


def plataforma(nombre):
    if nombre != "auto":
        return mm.Platform.getPlatformByName(nombre)
    for n in ("CUDA", "HIP", "OpenCL", "CPU"):
        try:
            return mm.Platform.getPlatformByName(n)
        except Exception:
            pass
    return mm.Platform.getPlatformByName("Reference")


ML = None   # configuración ML/MM (se rellena en main si se usa --ml)


def construir(prmtop, cutoff_nm):
    system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=cutoff_nm * u.nanometer,
                                 constraints=app.HBonds, rigidWater=True)
    if ML is not None:
        import potencial_ml
        system, _ = potencial_ml.convertir_a_mlmm(system, prmtop.topology, ML["posiciones"], ML["modelo"],
                                                  ML["indice"], ML["cuda_graphs"], ML["carpeta"],
                                                  verbose=ML.pop("verbose", False), doble=ML["doble"])
    return system


def añadir_restricciones(system, prmtop, posiciones, k_kcal_A2):
    f = mm.CustomExternalForce("k*periodicdistance(x,y,z,x0,y0,z0)^2")  # convenio AMBER: k·Δr²
    f.addGlobalParameter("k", k_kcal_A2 * KCAL_A2_TO_KJ_NM2)
    for p in ("x0", "y0", "z0"):
        f.addPerParticleParameter(p)
    n = 0
    for atom in prmtop.topology.atoms():
        if atom.residue.name in SOLVENTE:
            continue
        if atom.element is not None and atom.element.symbol == "H":
            continue
        x = posiciones[atom.index].value_in_unit(u.nanometer)
        f.addParticle(atom.index, list(x))
        n += 1
    system.addForce(f)
    return n


SOLVENTE = ("WAT", "HOH", "Na+", "Cl-", "NA", "CL", "K+")


def indices_soluto(topology):
    return [at.index for at in topology.atoms() if at.residue.name not in SOLVENTE]


def reporteros(sim, prefijo, cada_pasos_log, cada_pasos_dcd, total, dcd=True, subset=None):
    sim.reporters = []
    sim.reporters.append(app.StateDataReporter(prefijo + ".log", cada_pasos_log, step=True, time=True,
                         potentialEnergy=True, kineticEnergy=True, temperature=True, volume=True,
                         density=True, speed=True, separator="\t"))
    sim.reporters.append(app.StateDataReporter(sys.stdout, max(cada_pasos_log * 10, 1), step=True, time=True,
                         temperature=True, density=True, progress=True, remainingTime=True, speed=True,
                         totalSteps=total, separator="\t"))
    if dcd:
        if subset is None:
            sim.reporters.append(app.DCDReporter(prefijo + ".dcd", cada_pasos_dcd, enforcePeriodicBox=False))
        elif "atomSubset" in inspect.signature(app.DCDReporter.__init__).parameters:
            sim.reporters.append(app.DCDReporter(prefijo + ".dcd", cada_pasos_dcd, enforcePeriodicBox=False,
                                                 atomSubset=subset))
        else:
            # OpenMM antiguo (p. ej. 8.2, el que funciona con drivers CUDA 12.0-12.8): sin atomSubset
            sim.reporters.append(DCDSubconjunto(prefijo + ".dcd", cada_pasos_dcd, subset))


class DCDSubconjunto:
    """Escribe en DCD solo los átomos de `subset` (equivale a DCDReporter(atomSubset=...) de OpenMM >= 8.3)."""

    def __init__(self, archivo, cada, subset):
        self._cada = cada
        self._subset = list(subset)
        self._out = open(archivo, "wb")
        self._dcd = None

    def describeNextReport(self, simulation):
        pasos = self._cada - simulation.currentStep % self._cada
        return (pasos, True, False, False, False, False)  # posiciones sin envolver en la caja

    def report(self, simulation, state):
        if self._dcd is None:
            top = app.Topology()
            res = top.addResidue("SUB", top.addChain())
            for _ in self._subset:
                top.addAtom("X", None, res)
            top.setPeriodicBoxVectors(simulation.topology.getPeriodicBoxVectors())
            self._dcd = app.DCDFile(self._out, top, simulation.integrator.getStepSize(),
                                    simulation.currentStep, self._cada)
        pos = state.getPositions(asNumpy=True)[self._subset]
        self._dcd.writeModel(pos, periodicBoxVectors=state.getPeriodicBoxVectors())


def cerrar_reporteros(sim):
    for r in sim.reporters:
        f = getattr(r, "_out", None)
        if f is not None and f is not sys.stdout and hasattr(f, "close"):
            try:
                f.close()
            except Exception:
                pass
        dcd = getattr(r, "_dcd", None)
        if dcd is not None:
            r._dcd = None
    sim.reporters = []


def guardar(sim, prefijo):
    st = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=False)
    with open(prefijo + ".xml", "w") as fh:
        fh.write(mm.XmlSerializer.serialize(sim.context.getState(getPositions=True, getVelocities=True)))
    app.PDBFile.writeFile(sim.topology, st.getPositions(), open(prefijo + ".pdb", "w"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--sistema", default="sistema", help="carpeta con sistema.prmtop / sistema.inpcrd")
    ap.add_argument("--out", default="md")
    ap.add_argument("--etapas", default="min,cal,eq,prod")
    ap.add_argument("--temperatura", type=float, default=310.0)
    ap.add_argument("--presion-atm", type=float, default=1.0)
    ap.add_argument("--friccion", type=float, default=1.0, help="ps^-1")
    ap.add_argument("--dt", type=float, default=2.0, help="fs")
    ap.add_argument("--cutoff", type=float, default=11.0, help="Å")
    ap.add_argument("--k-restriccion", type=float, default=10.0, help="kcal/mol/Å² (minimización)")
    ap.add_argument("--min-iter", type=int, default=10000)
    ap.add_argument("--t-calentamiento", type=float, default=500.0, help="ps")
    ap.add_argument("--escalones", type=int, default=31, help="escalones de temperatura (10 K c/u)")
    ap.add_argument("--t-equilibrado", type=float, default=500.0, help="ps")
    ap.add_argument("--t-produccion", type=float, default=1_000_000.0, help="ps (1 µs en el artículo)")
    ap.add_argument("--segmentos", type=int, default=20,
                    help="la producción se divide en segmentos reanudables (20 x 50 ns en el artículo)")
    ap.add_argument("--max-segmentos", type=int, default=0,
                    help="correr como máximo N segmentos en esta ejecución (0 = todos)")
    ap.add_argument("--reintentos", type=int, default=2,
                    help="si un segmento produce NaN, repetirlo desde su estado inicial hasta N veces")
    ap.add_argument("--guardar-cada", type=float, default=10.0, help="ps entre fotogramas")
    ap.add_argument("--dcd-completo", action="store_true",
                    help="escribir agua e iones en la trayectoria (por defecto sólo proteína+ligando,\n"
                         "que se lee con sistema/complejo_seco.prmtop)")
    ap.add_argument("--plataforma", default="auto")
    ap.add_argument("--hilos", type=int, default=0)
    ap.add_argument("--semilla", type=int, default=2022)
    ap.add_argument("--ml", default="no", choices=("no", "ani2x", "ani1ccx"),
                    help="describir el ligando con un potencial de machine learning (ML/MM)")
    ap.add_argument("--ml-modelo-indice", type=int, default=-1,
                    help="-1 = conjunto de 8 redes (más exacto); 0-7 = una sola red (más rápido)")
    ap.add_argument("--ml-doble", action="store_true", help="evaluar la red en float64 (más lento)")
    ap.add_argument("--ml-cuda-graphs", action="store_true", help="usar CUDA graphs en la TorchForce")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    etapas = a.etapas.split(",")
    prmtop = app.AmberPrmtopFile(os.path.join(a.sistema, "sistema.prmtop"))
    inpcrd = app.AmberInpcrdFile(os.path.join(a.sistema, "sistema.inpcrd"))
    # una carpeta de trayectoria es de un solo método: no mezclar segmentos clásicos y ML/MM
    metodo = "clasico" if a.ml == "no" else f"mlmm-{a.ml}" + (f"-red{a.ml_modelo_indice}" if a.ml_modelo_indice >= 0 else "")
    marca = os.path.join(a.out, "metodo.txt")
    if os.path.exists(marca):
        previo = open(marca).read().strip()
        if previo != metodo:
            sys.exit(f"ERROR: {a.out} contiene una simulación '{previo}' y se pidió '{metodo}'. Usar otra carpeta --out.")
    elif any(n.endswith(".xml") for n in os.listdir(a.out)) and metodo != "clasico":
        sys.exit(f"ERROR: {a.out} ya tiene estados de una simulación clásica. Usar otra carpeta --out.")
    else:
        open(marca, "w").write(metodo + "\n")
    global ML
    if a.ml != "no":
        ML = {"modelo": a.ml, "indice": a.ml_modelo_indice, "cuda_graphs": a.ml_cuda_graphs, "doble": a.ml_doble,
              "posiciones": inpcrd.positions, "carpeta": a.out, "verbose": True}
        print(f"Método: ML/MM, ligando con {a.ml.upper()} (proteína, agua e iones con ff19SB/TIP3P)")
    T = a.temperatura * u.kelvin
    dt = a.dt * u.femtoseconds
    ps2steps = lambda t_ps: int(round(t_ps / (a.dt / 1000.0)))
    plat = plataforma(a.plataforma)
    props = {}
    if plat.getName() == "CPU" and a.hilos > 0:
        props["Threads"] = str(a.hilos)
    if plat.getName() in ("CUDA", "OpenCL", "HIP"):
        props["Precision"] = "mixed"
    print(f"Plataforma: {plat.getName()} {props}")
    print(f"Átomos: {prmtop.topology.getNumAtoms()}  |  caja: {inpcrd.boxVectors}")

    def nueva_sim(system, temp=T):
        integ = mm.LangevinMiddleIntegrator(temp, a.friccion / u.picosecond, dt)
        integ.setRandomNumberSeed(a.semilla)
        integ.setConstraintTolerance(1e-6)
        return app.Simulation(prmtop.topology, system, integ, plat, props)

    out = lambda n: os.path.join(a.out, n)
    subset = None if a.dcd_completo else indices_soluto(prmtop.topology)

    # -------------------------------------------------------------- minimización
    if "min" in etapas:
        print("\n=== Minimización (restricción %.1f kcal/mol/Å², máx. %d iteraciones) ===" % (a.k_restriccion, a.min_iter))
        system = construir(prmtop, a.cutoff / 10.0)
        n = añadir_restricciones(system, prmtop, inpcrd.positions, a.k_restriccion)
        print(f"  átomos restringidos: {n}")
        sim = nueva_sim(system)
        sim.context.setPositions(inpcrd.positions)
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
        t0 = time.time()
        sim.minimizeEnergy(tolerance=10 * u.kilojoule_per_mole / u.nanometer, maxIterations=a.min_iter)
        e1 = sim.context.getState(getEnergy=True).getPotentialEnergy()
        print(f"  E inicial = {e0}\n  E final   = {e1}   ({time.time()-t0:.0f} s)")
        with open(out("01_min.txt"), "w") as fh:
            fh.write(f"E_inicial {e0}\nE_final {e1}\n")
        guardar(sim, out("01_min"))

    # -------------------------------------------------------------- calentamiento
    if "cal" in etapas:
        print(f"\n=== Calentamiento 0 → {a.temperatura:.0f} K en {a.t_calentamiento} ps (NVT, sin restricciones) ===")
        system = construir(prmtop, a.cutoff / 10.0)
        sim = nueva_sim(system, 0.0 * u.kelvin)
        sim.loadState(out("01_min.xml"))
        total = ps2steps(a.t_calentamiento)
        por_escalon = max(total // a.escalones, 1)
        reporteros(sim, out("02_cal"), max(ps2steps(1.0), 1), ps2steps(a.guardar_cada), total, subset=subset)
        sim.context.setVelocitiesToTemperature(1.0 * u.kelvin, a.semilla)
        for i in range(1, a.escalones + 1):
            Ti = a.temperatura * i / a.escalones
            sim.integrator.setTemperature(Ti * u.kelvin)
            sim.step(por_escalon if i < a.escalones else total - por_escalon * (a.escalones - 1))
        guardar(sim, out("02_cal"))

    # -------------------------------------------------------------- equilibrado
    barostato = lambda s: s.addForce(mm.MonteCarloBarostat(a.presion_atm * u.atmosphere, T, 25))
    if "eq" in etapas:
        print(f"\n=== Equilibrado {a.t_equilibrado} ps (NPT, {a.temperatura:.0f} K, {a.presion_atm} atm) ===")
        system = construir(prmtop, a.cutoff / 10.0); barostato(system)
        sim = nueva_sim(system)
        sim.loadState(out("02_cal.xml"))
        sim.integrator.setTemperature(T)
        total = ps2steps(a.t_equilibrado)
        reporteros(sim, out("03_eq"), max(ps2steps(1.0), 1), ps2steps(a.guardar_cada), total, subset=subset)
        sim.step(total)
        guardar(sim, out("03_eq"))

    # -------------------------------------------------------------- producción
    if "prod" in etapas:
        nseg = max(a.segmentos, 1)
        t_seg = a.t_produccion / nseg
        print(f"\n=== Producción {a.t_produccion} ps en {nseg} segmento(s) de {t_seg} ps "
              f"(NPT, {a.temperatura:.0f} K, {a.presion_atm} atm) ===")
        system = construir(prmtop, a.cutoff / 10.0); barostato(system)
        sim = nueva_sim(system)
        hechos = 0
        for i in range(1, nseg + 1):
            if a.max_segmentos and hechos >= a.max_segmentos:
                print(f"  (se alcanzó --max-segmentos={a.max_segmentos}; volver a ejecutar para continuar)")
                break
            pref = out(f"04_prod_{i:03d}")
            if os.path.exists(pref + ".xml"):
                print(f"  segmento {i} ya terminado, se omite")
                continue
            prev = out("03_eq.xml") if i == 1 else out(f"04_prod_{i-1:03d}.xml")
            total = ps2steps(t_seg)
            for intento in range(1, a.reintentos + 2):
                sim.loadState(prev)
                sim.integrator.setTemperature(T)
                sim.currentStep = ps2steps(t_seg) * (i - 1)
                sim.context.setTime(t_seg * (i - 1) * u.picosecond)
                print(f"  -- segmento {i}/{nseg}  ({t_seg*(i-1):.1f} → {t_seg*i:.1f} ps)"
                      + (f"  [reintento {intento-1}]" if intento > 1 else ""))
                reporteros(sim, pref, ps2steps(a.guardar_cada), ps2steps(a.guardar_cada), total, subset=subset)
                sim.reporters.append(app.CheckpointReporter(pref + ".chk", max(ps2steps(min(100.0, t_seg)), 1)))
                try:
                    sim.step(total)
                    break
                except mm.OpenMMException as e:
                    cerrar_reporteros(sim)
                    if "NaN" not in str(e) or intento > a.reintentos:
                        raise
                    # el ruido de Langevin sigue avanzando, así que repetir desde el mismo
                    # estado da una trayectoria distinta; queda registrado en el log
                    msg = f"segmento {i}: NaN en t≈{sim.currentStep*a.dt/1000:.1f} ps; se repite desde {os.path.basename(prev)}"
                    print("  !! " + msg)
                    with open(out("incidencias.log"), "a") as fh:
                        fh.write(msg + "\n")
            cerrar_reporteros(sim)
            guardar(sim, pref)
            hechos += 1

    print("\nTerminado.")


if __name__ == "__main__":
    main()
