#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
potencial_ml.py — convierte el sistema clásico (ff19SB + GAFF) en un sistema ML/MM: el ligando
(residuo LIG) se describe con el potencial de red neuronal ANI-2x (TorchANI) y la proteína, el agua
y los iones siguen con el campo de fuerza.

Acoplamiento mecánico (mechanical embedding), el mismo esquema que usa openmm-ml:
  - Se anulan los enlaces, ángulos y torsiones (incluidos los impropios) con todos sus átomos en el
    ligando. Se dejan con constante 0 para que OpenMM siga viendo el ligando como una sola molécula.
  - Se anulan las interacciones no enlazadas ligando–ligando (excepciones con carga y ε = 0).
  - Las interacciones ligando–entorno (Lennard-Jones y cargas GAFF/AM1-BCC) siguen siendo clásicas.
  - Se añade una TorchForce con la energía interna del ligando calculada por ANI-2x.
Se mantienen las restricciones de los enlaces con H del ligando, como en la versión clásica, para
conservar el paso de 2 fs.

Limitaciones de ANI-2x: solo moléculas neutras con H, C, N, O, F, S y Cl. El script se detiene si
el ligando tiene carga (por ejemplo, Pallidine protonada).

Verificación (no simula nada, solo construye el sistema y evalúa la energía y las fuerzas):
  python potencial_ml.py --sistema ../sistema
"""
import argparse, json, os, sys
import numpy as np
import openmm as mm
from openmm import app, unit as u

ELEMENTOS_ANI2X = {"H", "C", "N", "O", "F", "S", "Cl"}
MODELOS = ("ani2x", "ani1ccx")


def atomos_ligando(topology, resname="LIG"):
    return [a.index for a in topology.atoms() if a.residue.name == resname]


def carga_ligando(system, idx):
    nb = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
    return sum(nb.getParticleParameters(i)[0].value_in_unit(u.elementary_charge) for i in idx)


def verificar_ligando(system, topology, idx, modelo):
    if not idx:
        sys.exit("ERROR: no hay ningún residuo LIG en la topología.")
    atomos = list(topology.atoms())
    elem = sorted({atomos[i].element.symbol for i in idx})
    raros = [e for e in elem if e not in ELEMENTOS_ANI2X or (modelo == "ani1ccx" and e not in "HCNO")]
    if raros:
        sys.exit(f"ERROR: {modelo} no admite los elementos {raros} del ligando.")
    q = carga_ligando(system, idx)
    if abs(q) > 0.5:
        sys.exit(f"ERROR: el ligando tiene carga {q:+.2f}. {modelo.upper()} solo está entrenado con moléculas "
                 "neutras, así que no puede describir esta especie (por ejemplo, Pallidine protonada). "
                 "Usar la forma neutra o la dinámica clásica.")
    return elem, q


try:
    import torch
    import torchani
except ImportError as e:  # el resto del repositorio no necesita PyTorch
    sys.exit(f"ERROR: la versión ML necesita PyTorch y TorchANI (entorno md-braf-ml): {e}")


class EnergiaANI(torch.nn.Module):
    """Módulo TorchScript: posiciones (nm) y caja (nm) -> energía interna del ligando (kJ/mol)."""

    def __init__(self, red, especies, indices, doble: bool = False):
        super().__init__()
        self.red = red
        self.doble = doble
        self.register_buffer("especies", especies)
        self.register_buffer("indices", indices)
        self.register_buffer("e_ref", torch.zeros(1, dtype=torch.float64))
        self.escala = float(torchani.units.hartree2kjoulemol(1.0))

    def forward(self, positions: torch.Tensor, boxvectors: torch.Tensor) -> torch.Tensor:
        if self.doble:
            x = positions.index_select(0, self.indices).to(torch.float64) * 10.0   # nm -> Å
            b = boxvectors.to(torch.float64) * 10.0
        else:
            x = positions.index_select(0, self.indices).to(torch.float32) * 10.0
            b = boxvectors.to(torch.float32) * 10.0
        # imagen mínima respecto al primer átomo: el ligando nunca queda partido por la caja
        d = x - x[0:1]
        d = d - torch.round(d @ torch.linalg.inv(b)) @ b
        x = x[0:1] + d
        e = self.red((self.especies, x.unsqueeze(0))).energies.to(torch.float64) * self.escala
        return (e - self.e_ref).sum()


def _modulo_ani(modelo, especies, idx, modelo_indice, doble=False):
    kw = {"periodic_table_index": True}
    if modelo_indice is not None and modelo_indice >= 0:
        kw["model_index"] = modelo_indice
    red = getattr(torchani.models, {"ani2x": "ANI2x", "ani1ccx": "ANI1ccx"}[modelo])(**kw)
    if doble:
        red = red.to(torch.float64)
    return EnergiaANI(red, torch.tensor([especies], dtype=torch.int64), torch.tensor(idx, dtype=torch.int64), doble)


def convertir_a_mlmm(system, topology, posiciones, modelo="ani2x", modelo_indice=None,
                     cuda_graphs=False, carpeta=".", verbose=True, doble=False):
    """Devuelve una copia de `system` con el ligando descrito por `modelo` (ML/MM)."""
    from openmmtorch import TorchForce

    idx = atomos_ligando(topology)
    elem, q = verificar_ligando(system, topology, idx, modelo)
    s = set(idx)
    new = mm.XmlSerializer.deserialize(mm.XmlSerializer.serialize(system))
    anulados = {"enlaces": 0, "ángulos": 0, "torsiones": 0, "pares no enlazados": 0}
    for f in new.getForces():
        if isinstance(f, mm.HarmonicBondForce):
            for i in range(f.getNumBonds()):
                p1, p2, r0, k = f.getBondParameters(i)
                if p1 in s and p2 in s:
                    f.setBondParameters(i, p1, p2, r0, 0.0); anulados["enlaces"] += 1
        elif isinstance(f, mm.HarmonicAngleForce):
            for i in range(f.getNumAngles()):
                p1, p2, p3, th, k = f.getAngleParameters(i)
                if p1 in s and p2 in s and p3 in s:
                    f.setAngleParameters(i, p1, p2, p3, th, 0.0); anulados["ángulos"] += 1
        elif isinstance(f, mm.PeriodicTorsionForce):
            for i in range(f.getNumTorsions()):
                p1, p2, p3, p4, n, ph, k = f.getTorsionParameters(i)
                if p1 in s and p2 in s and p3 in s and p4 in s:
                    f.setTorsionParameters(i, p1, p2, p3, p4, n, ph, 0.0); anulados["torsiones"] += 1
        elif isinstance(f, mm.NonbondedForce):
            for a in range(len(idx)):
                for b in range(a):
                    f.addException(idx[a], idx[b], 0.0, 1.0, 0.0, True)
                    anulados["pares no enlazados"] += 1
        elif isinstance(f, mm.CustomNonbondedForce):
            excl = {tuple(sorted(f.getExclusionParticles(i))) for i in range(f.getNumExclusions())}
            for a in range(len(idx)):
                for b in range(a):
                    if tuple(sorted((idx[a], idx[b]))) not in excl:
                        f.addExclusion(idx[a], idx[b])
        elif isinstance(f, mm.CMAPTorsionForce):
            pass  # solo proteína
        elif isinstance(f, (mm.CMMotionRemover, mm.MonteCarloBarostat, mm.CustomExternalForce)):
            pass
        else:
            sys.exit(f"ERROR: fuerza no prevista en el sistema: {f.__class__.__name__}")

    atomos = list(topology.atoms())
    especies = [atomos[i].element.atomic_number for i in idx]
    mod = _modulo_ani(modelo, especies, idx, modelo_indice, doble)
    # referencia de energía: ANI en la geometría inicial (solo desplaza la energía potencial
    # registrada para que quede en la misma escala que la clásica; no cambia las fuerzas)
    with torch.no_grad():
        x = torch.tensor(np.array(posiciones.value_in_unit(u.nanometer), dtype=np.float64))
        bv = np.array([[v.value_in_unit(u.nanometer) for v in vec] for vec in topology.getPeriodicBoxVectors()]) \
            if topology.getPeriodicBoxVectors() is not None else np.eye(3) * 100.0
        mod.e_ref.copy_(mod(x, torch.tensor(bv)).reshape(1))
    os.makedirs(carpeta, exist_ok=True)
    ruta = os.path.join(carpeta, f"{modelo}_ligando.pt")
    torch.jit.script(mod).save(ruta)
    tf = TorchForce(ruta)
    tf.setUsesPeriodicBoundaryConditions(True)
    if cuda_graphs:
        try:
            tf.setProperty("useCUDAGraphs", "true")
        except Exception:
            print("  (esta versión de openmm-torch no admite CUDA graphs; se ignora)")
    new.addForce(tf)
    info = {"modelo": modelo, "modelo_indice": modelo_indice, "atomos_ligando": len(idx), "elementos": elem,
            "carga_ligando": round(q, 3), "terminos_anulados": anulados, "restricciones_H_ligando": "se mantienen",
            "acoplamiento": "mecánico (ligando–entorno con GAFF/AM1-BCC + LJ)",
            "precision_red": "float64" if doble else "float32", "e_ref_kJ_mol": float(mod.e_ref.item())}
    if verbose:
        print(f"  ML/MM: {modelo}{' (red %d)' % modelo_indice if modelo_indice is not None and modelo_indice >= 0 else ' (conjunto de 8 redes)'}"
              f" para {len(idx)} átomos del ligando ({', '.join(elem)}; carga {q:+.2f})")
        print(f"         términos clásicos anulados: {anulados}")
    return new, info


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--sistema", default="../sistema")
    ap.add_argument("--ml", default="ani2x", choices=MODELOS)
    ap.add_argument("--ml-modelo-indice", type=int, default=-1, help="-1 = conjunto de 8 redes")
    ap.add_argument("--plataforma", default="auto")
    ap.add_argument("--out", default=None, help="carpeta para guardar el modelo TorchScript (por defecto, temporal)")
    a = ap.parse_args()
    import tempfile
    prmtop = app.AmberPrmtopFile(os.path.join(a.sistema, "sistema.prmtop"))
    inpcrd = app.AmberInpcrdFile(os.path.join(a.sistema, "sistema.inpcrd"))
    clas = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.1 * u.nanometer,
                               constraints=app.HBonds, rigidWater=True)
    carpeta = a.out or tempfile.mkdtemp()
    for i, f in enumerate(clas.getForces()):
        f.setForceGroup(i % 31)
    mlsys, info = convertir_a_mlmm(clas, prmtop.topology, inpcrd.positions, a.ml, a.ml_modelo_indice, carpeta=carpeta)
    mlsys.getForce(mlsys.getNumForces() - 1).setForceGroup(31)
    if a.plataforma == "auto":
        for n in ("CUDA", "HIP", "OpenCL", "CPU"):
            try:
                plat = mm.Platform.getPlatformByName(n); break
            except Exception:
                pass
    else:
        plat = mm.Platform.getPlatformByName(a.plataforma)
    print(f"  Plataforma OpenMM: {plat.getName()}  |  PyTorch {torch.__version__}, CUDA disponible: {torch.cuda.is_available()}")
    ctx = mm.Context(mlsys, mm.VerletIntegrator(0.001), plat)
    ctx.setPositions(inpcrd.positions); ctx.setPeriodicBoxVectors(*inpcrd.boxVectors)
    st = ctx.getState(getEnergy=True, getForces=True, groups={31})
    fl = st.getForces(asNumpy=True).value_in_unit(u.kilojoule_per_mole / u.nanometer)[atomos_ligando(prmtop.topology)]
    print(f"  Energía ANI (relativa a la geometría inicial): {st.getPotentialEnergy()}")
    print(f"  |F| máxima sobre el ligando: {np.abs(fl).max():.1f} kJ/mol/nm  (finitas: {np.isfinite(fl).all()})")
    ctx2 = mm.Context(clas, mm.VerletIntegrator(0.001), plat)
    ctx2.setPositions(inpcrd.positions); ctx2.setPeriodicBoxVectors(*inpcrd.boxVectors)
    e_cl = ctx2.getState(getEnergy=True).getPotentialEnergy()
    e_ml = ctx.getState(getEnergy=True).getPotentialEnergy()
    print(f"  Energía total: clásica {e_cl}  |  ML/MM {e_ml}")
    print(json.dumps(info, ensure_ascii=False, indent=1))
    print("Verificación correcta: el sistema ML/MM se puede simular.")


if __name__ == "__main__":
    main()
