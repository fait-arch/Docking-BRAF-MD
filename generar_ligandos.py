#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Genera los PDB de los ligandos en su estado de protonación a pH 7.4.

    python3 generar_ligandos.py

Los PDB que venían de las bases de datos están dibujados en forma neutra. A
pH 7.4 varios de estos compuestos no lo están, y la carga no es un número que
se pueda cambiar en el YAML: `mol_ligando()` deduce los enlaces de la geometría
contando los hidrógenos del archivo, así que el PDB tiene que llevar el protón
puesto (o quitado).

Cada SMILES de partida es el que el propio pipeline dedujo del PDB original
(simulaciones/<X>/sistema/sistema_info.json), con un solo átomo cambiado. Así
se conserva la estereoquímica y no se introduce otra molécula por error.

Cada archivo se verifica releyéndolo con el mismo mol_ligando() del pipeline y
comparando el InChIKey contra el del SMILES de referencia.
"""
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")
RAIZ = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(RAIZ, "Data")

# (archivo, SMILES, carga, nota sobre el porqué)
LIGANDOS = [
    ("Quinine_72243208_cation_ligand.pdb",
     "C=C[C@H]1C[NH+]2CC[C@H]1C[C@H]2[C@H](O)c1ccnc2ccc(OC)cc12", 1,
     "N de la quinuclidina protonado. pKa 8.4 (Roy et al., Br J Clin "
     "Pharmacol 2002;54:604-609) -> 91 % cationica a pH 7.4. El N "
     "quinolinico, aromatico, se deja neutro."),

    ("NotholaenicAcid_78213884_anion_ligand.pdb",
     "COc1ccc(CCc2cc(OC)cc(O)c2C(=O)[O-])cc1", -1,
     "Carboxilo desprotonado. El grupo esta confirmado (Cioffi et al., "
     "Molecules 2011;16:2527-2541: 5-hidroxi-3,12-dimetoxi-6-carboxibibencilo) "
     "pero NO hay pKa experimental publicado: este estado es una inferencia "
     "del andamiaje, por eso se conserva tambien la forma neutra."),

    ("Pallidine_24543347_cation_ligand.pdb",
     "COC1=C[C@@]23CC[NH+](C)[C@@H](Cc4cc(O)c(OC)cc42)C3=CC1=O", 1,
     "Amina terciaria protonada. Es el SMILES que ya estaba declarado en "
     "config/compuestos.yaml para pallidine_protonada, que hasta ahora no "
     "tenia PDB propio y por eso fallaba."),
]


def como_el_pipeline(pdb, carga):
    """Copia literal de mol_ligando() de 01_preparar_sistema.py."""
    m = Chem.MolFromPDBFile(pdb, removeHs=False, sanitize=False,
                            proximityBonding=True)
    rdDetermineBonds.DetermineBondOrders(m, charge=carga)
    Chem.SanitizeMol(m)
    return m


def generar(nombre, smiles, carga):
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        return f"SMILES invalido"
    if Chem.GetFormalCharge(ref) != carga:
        return f"el SMILES tiene carga {Chem.GetFormalCharge(ref)}, no {carga}"
    clave_ref = Chem.MolToInchiKey(ref)

    m = Chem.AddHs(ref)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 2022                      # misma semilla que el pipeline
    if AllChem.EmbedMolecule(m, ps) != 0:
        return "no se pudo generar la conformacion 3D"
    if AllChem.MMFFHasAllMoleculeParams(m):
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    else:
        AllChem.UFFOptimizeMolecule(m, maxIters=2000)

    cont = {}
    for a in m.GetAtoms():
        s = a.GetSymbol()
        cont[s] = cont.get(s, 0) + 1
        ri = Chem.AtomPDBResidueInfo()
        ri.SetResidueName("LIG")
        ri.SetResidueNumber(1)
        ri.SetChainId("A")
        ri.SetName(f"{s}{cont[s]}".ljust(4)[:4])
        ri.SetIsHeteroAtom(True)
        a.SetMonomerInfo(ri)

    ruta = os.path.join(DATA, nombre)
    Chem.MolToPDBFile(m, ruta, flavor=4)

    # -- la verificacion que importa: releerlo como lo hara el pipeline -----
    try:
        leido = como_el_pipeline(ruta, carga)
    except Exception as e:
        return f"el pipeline NO puede leerlo: {type(e).__name__}: {e}"
    clave = Chem.MolToInchiKey(Chem.RemoveHs(leido))
    if clave != clave_ref:
        return f"no coincide al releer: {clave} != {clave_ref}"
    return (CalcMolFormula(leido), Chem.GetFormalCharge(leido), clave)


def main():
    fallos = 0
    for nombre, smiles, carga, nota in LIGANDOS:
        r = generar(nombre, smiles, carga)
        if isinstance(r, str):
            print(f"  FALLA  {nombre}\n         {r}")
            fallos += 1
            continue
        formula, q, clave = r
        print(f"  OK     {nombre}")
        print(f"         {formula}   carga {q:+d}   {clave}")
        print(f"         {nota}")
        print()
    sys.exit(1 if fallos else 0)


if __name__ == "__main__":
    main()
