#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
01_preparar_sistema.py  —  Macaridine (PubChem 94717379) + BRAF V600E

Prepara el complejo proteína-ligando siguiendo la sección 2.1.3 de
RSC Adv. 2022, 12, 30181 (d2ra05751d):

  * Proteína: AMBER ff19SB, protonación a pH fisiológico (PDB2PQR + PROPKA, pH 7.4)
  * Ligando:  GAFF + cargas AM1-BCC (antechamber / parmchk2)
  * Pose inicial: docking con AutoDock Vina (caja 20x20x20 Å, exhaustiveness = 100,
    num_modes = 100), como en la sección 2.1.2 del artículo
  * Solvatación: agua TIP3P en caja octaédrica truncada, 12 Å al borde
  * Iones: neutralización + NaCl 0.15 M según SLTCAP (Machado & Pantano)

Uso:
  python 01_preparar_sistema.py --receptor ../entrada/receptor.pdb \
         --ligando ../entrada/ligand.pdb --out ../sistema

Requiere: ambertools (tleap, antechamber, parmchk2), pdb2pqr, vina, openbabel, rdkit.
"""
import argparse, json, math, os, re, shutil, subprocess, sys
import numpy as np

# Residuos que tapizan el bolsillo ATP / sitio de dabrafenib en BRAF
# (lazo P, Lys483, gatekeeper Thr529, bisagra Trp531-Cys532, DFG Asp594)
RES_BOLSILLO = [463, 471, 481, 483, 505, 514, 529, 531, 532, 583, 594]


def sh(cmd, log=None, cwd=None):
    print(f"  $ {cmd}")
    r = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    if log:
        with open(log, "a") as fh:
            fh.write(f"$ {cmd}\n{r.stdout}\n{r.stderr}\n")
    if r.returncode != 0:
        print(r.stdout[-3000:]); print(r.stderr[-3000:])
        raise SystemExit(f"ERROR ejecutando: {cmd}")
    return r.stdout


def leer_pdb_atoms(path):
    out = []
    for l in open(path):
        if l.startswith(("ATOM", "HETATM")):
            out.append(l)
    return out


def recortar_receptor(pdb_in, pdb_out, ini, fin):
    n = 0
    with open(pdb_out, "w") as fh:
        for l in leer_pdb_atoms(pdb_in):
            r = int(l[22:26])
            if ini <= r <= fin and l.startswith("ATOM"):
                fh.write(l); n += 1
        fh.write("TER\nEND\n")
    return n


def limpiar_pdb2pqr(pdb_in, pdb_out):
    """Quita hidrógenos (tleap los reconstruye) y normaliza nombres de residuo
    terminales que PDB2PQR escribe con prefijo N/C (p. ej. NSER, CVAL)."""
    std = {"ALA","ARG","ASN","ASP","ASH","CYS","CYX","CYM","GLN","GLU","GLH","GLY",
           "HID","HIE","HIP","HIS","ILE","LEU","LYS","LYN","MET","PHE","PRO","SER",
           "THR","TRP","TYR","VAL"}
    prot = {}
    with open(pdb_out, "w") as fh:
        for l in open(pdb_in):
            if not l.startswith("ATOM"):
                continue
            name = l[12:16].strip()
            elem = l[76:78].strip()
            if elem == "H" or (not elem and name.startswith("H")) or re.match(r"^\d*H", name):
                continue
            resn = l[17:21].strip()
            if resn not in std and len(resn) == 4 and resn[1:] in std | {"HIS"}:
                resn = resn[1:]
            if name == "OXT":
                pass
            prot[int(l[22:26])] = resn
            l = l[:17] + f"{resn:>3s} " + l[21:]
            fh.write(l)
        fh.write("TER\nEND\n")
    return prot


def centro_bolsillo(pdb, residuos):
    xyz = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
           for l in leer_pdb_atoms(pdb) if int(l[22:26]) in residuos]
    return np.mean(np.array(xyz), axis=0)


def mol_ligando(pdb_lig, carga):
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    m = Chem.MolFromPDBFile(pdb_lig, removeHs=False, sanitize=False, proximityBonding=True)
    rdDetermineBonds.DetermineBondOrders(m, charge=carga)
    Chem.SanitizeMol(m)
    return m


def sltcap(n_agua, q, conc=0.15):
    """Número de cationes/aniones según SLTCAP (Schmit et al. 2018, JCTC 14, 1823).
    N0 = n_agua * c / 55.5 M ;  N± = N0*(sqrt(1+y^2) ∓ y), y = Q/(2 N0)."""
    n0 = n_agua * conc / 55.5
    y = q / (2.0 * n0)
    npos = n0 * (math.sqrt(1 + y * y) - y)
    nneg = n0 * (math.sqrt(1 + y * y) + y)
    npos, nneg = int(round(npos)), int(round(nneg))
    # garantizar neutralidad exacta
    if npos - nneg + q != 0:
        nneg = npos + int(round(q))
    return npos, nneg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--ligando", required=True)
    ap.add_argument("--out", default="sistema")
    ap.add_argument("--residuos", default="448-723",
                    help="Rango de residuos de la proteína a conservar. El artículo usa el dominio\n"
                         "quinasa de BRAF V600E (PDB 4XV2). Use '1-766' para la proteína completa.")
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--carga-ligando", type=int, default=0)
    ap.add_argument("--sin-docking", action="store_true",
                    help="Usar las coordenadas del ligando tal como vienen (sin Vina).")
    ap.add_argument("--exhaustiveness", type=int, default=100)
    ap.add_argument("--caja", type=float, default=12.0, help="Distancia soluto-borde (Å)")
    ap.add_argument("--sal", type=float, default=0.15, help="NaCl (M)")
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    out = os.path.abspath(a.out)
    log = os.path.join(out, "preparacion.log")
    open(log, "w").close()
    info = {}
    ini, fin = [int(x) for x in a.residuos.split("-")]

    # ---------------------------------------------------------------- receptor
    print("\n[1] Receptor: recorte y protonación (PDB2PQR/PROPKA)")
    rec_trim = os.path.join(out, "receptor_recortado.pdb")
    n = recortar_receptor(a.receptor, rec_trim, ini, fin)
    print(f"    residuos {ini}-{fin}: {n} átomos pesados")
    rec_pqr = os.path.join(out, "receptor.pqr")
    rec_h = os.path.join(out, "receptor_H.pdb")
    sh(f"pdb2pqr30 --ff=AMBER --ffout=AMBER --titration-state-method=propka "
       f"--with-ph={a.ph} --keep-chain --pdb-output {rec_h} {rec_trim} {rec_pqr}", log)
    rec_noh = os.path.join(out, "receptor_amber.pdb")
    nombres = limpiar_pdb2pqr(rec_h, rec_noh)
    estados = {k: v for k, v in nombres.items() if v in ("HID","HIE","HIP","ASH","GLH","LYN","CYX","CYM")}
    info["protonacion"] = {str(k): v for k, v in sorted(estados.items())}
    print(f"    estados no estándar: {info['protonacion']}")

    # ---------------------------------------------------------------- ligando
    print("\n[2] Ligando: perceptión de enlaces")
    from rdkit import Chem
    from rdkit.Chem import AllChem
    lig = mol_ligando(a.ligando, a.carga_ligando)
    info["ligando_smiles"] = Chem.MolToSmiles(Chem.RemoveHs(lig))
    print(f"    SMILES: {info['ligando_smiles']}")

    centro = centro_bolsillo(rec_noh, RES_BOLSILLO)
    lig_xyz = lig.GetConformer().GetPositions()
    info["distancia_ligando_original_al_bolsillo_A"] = float(np.linalg.norm(lig_xyz.mean(0) - centro))

    if a.sin_docking:
        lig_pose = lig
    else:
        print("\n[3] Docking AutoDock Vina (caja 20 Å, exhaustiveness=%d)" % a.exhaustiveness)
        rec_pdbqt = os.path.join(out, "receptor.pdbqt")
        sh(f"obabel {rec_h} -xr -O {rec_pdbqt}", log)
        # ligando de partida: geometría propia (RDKit), lejos de la proteína
        lig_in = Chem.Mol(lig)
        AllChem.MMFFOptimizeMolecule(lig_in)
        sdf_in = os.path.join(out, "ligando_entrada.sdf")
        Chem.MolToMolFile(lig_in, sdf_in)
        lig_pdbqt = os.path.join(out, "ligando.pdbqt")
        sh(f"obabel {sdf_in} -O {lig_pdbqt}", log)
        from vina import Vina
        v = Vina(sf_name="vina", cpu=0, seed=42)
        v.set_receptor(rec_pdbqt)
        v.set_ligand_from_file(lig_pdbqt)
        v.compute_vina_maps(center=list(map(float, centro)), box_size=[20, 20, 20])
        v.dock(exhaustiveness=a.exhaustiveness, n_poses=100)
        poses_out = os.path.join(out, "vina_poses.pdbqt")
        v.write_poses(poses_out, n_poses=100, overwrite=True)
        e = v.energies(n_poses=100)
        info["vina_mejor_kcal_mol"] = float(e[0][0])
        info["vina_n_poses"] = int(len(e))
        print(f"    mejor pose: {e[0][0]:.2f} kcal/mol ({len(e)} poses)")
        # mejor pose -> molécula completa con H
        best = os.path.join(out, "vina_pose1.pdbqt")
        v.write_poses(best, n_poses=1, overwrite=True)
        best_sdf = os.path.join(out, "pose1.sdf")
        sh(f"obabel {best} -O {best_sdf}", log)
        pose = Chem.MolFromMolFile(best_sdf, removeHs=False, sanitize=False)
        tmpl = Chem.RemoveHs(lig)
        pose_heavy = Chem.RemoveHs(pose, sanitize=False)
        pose_heavy = AllChem.AssignBondOrdersFromTemplate(tmpl, pose_heavy)
        lig_pose = Chem.AddHs(pose_heavy, addCoords=True)
        Chem.MolToMolFile(lig_pose, os.path.join(out, "ligando_docking.sdf"))

    # ---------------------------------------------------------------- GAFF
    print("\n[4] Parámetros GAFF + AM1-BCC (antechamber)")
    lig_sdf = os.path.join(out, "ligando_pose.sdf")
    Chem.MolToMolFile(lig_pose, lig_sdf)
    sh(f"antechamber -i ligando_pose.sdf -fi sdf -o LIG.mol2 -fo mol2 -c bcc "
       f"-nc {a.carga_ligando} -at gaff -rn LIG -pf y", log, cwd=out)
    sh("parmchk2 -i LIG.mol2 -f mol2 -o LIG.frcmod -s gaff", log, cwd=out)

    # distancia mínima ligando-proteína de la pose final
    pxyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in leer_pdb_atoms(rec_noh)])
    lxyz = lig_pose.GetConformer().GetPositions()
    dmin = float(np.min(np.linalg.norm(lxyz[:, None, :] - pxyz[None], axis=2)))
    info["distancia_minima_ligando_proteina_A"] = round(dmin, 2)
    print(f"    distancia mínima ligando-proteína (pose final): {dmin:.2f} Å")

    # ---------------------------------------------------------------- tleap
    print("\n[5] tleap: ff19SB + GAFF + TIP3P, octaedro truncado, SLTCAP")
    base = (
        "source leaprc.protein.ff19SB\n"
        "source leaprc.gaff\n"
        "source leaprc.water.tip3p\n"
        "loadamberparams LIG.frcmod\n"
        "LIG = loadmol2 LIG.mol2\n"
        "prot = loadpdb receptor_amber.pdb\n"
        "complex = combine {prot LIG}\n"
    )
    with open(os.path.join(out, "tleap_seco.in"), "w") as fh:
        fh.write(base + "check complex\ncharge complex\nsavepdb complex complejo_seco.pdb\n"
                 "saveamberparm complex complejo_seco.prmtop complejo_seco.inpcrd\n"
                 f"solvateOct complex TIP3PBOX {a.caja}\nquit\n")
    o = sh("tleap -f tleap_seco.in", log, cwd=out)
    q = float(re.findall(r"Total unperturbed charge:\s+(-?[\d.]+)", o)[-1])
    nwat = int(re.findall(r"Added (\d+) residues", o)[-1])
    npos, nneg = sltcap(nwat, round(q), a.sal)
    info.update({"carga_soluto": q, "aguas_iniciales": nwat, "Na+": npos, "Cl-": nneg})
    print(f"    carga del soluto = {q:+.0f}, aguas = {nwat}  ->  Na+ = {npos}, Cl- = {nneg}")

    with open(os.path.join(out, "tleap.in"), "w") as fh:
        fh.write(base +
                 f"solvateOct complex TIP3PBOX {a.caja}\n"
                 f"addIonsRand complex Na+ {npos} Cl- {nneg} 5.0\n"
                 "charge complex\n"
                 "saveamberparm complex sistema.prmtop sistema.inpcrd\n"
                 "savepdb complex sistema.pdb\nquit\n")
    o = sh("tleap -f tleap.in", log, cwd=out)
    qf = float(re.findall(r"Total unperturbed charge:\s+(-?[\d.]+)", o)[-1])

    import parmed as pmd
    s = pmd.load_file(os.path.join(out, "sistema.prmtop"), os.path.join(out, "sistema.inpcrd"))
    info.update({"carga_final": qf, "n_atomos": len(s.atoms),
                 "aguas_finales": sum(1 for r in s.residues if r.name == "WAT"),
                 "caja_A_grados": [round(x, 3) for x in s.box]})
    with open(os.path.join(out, "sistema_info.json"), "w") as fh:
        json.dump(info, fh, indent=2, ensure_ascii=False)
    print("\n" + json.dumps(info, indent=2, ensure_ascii=False))
    print("\nListo: sistema.prmtop / sistema.inpcrd")


if __name__ == "__main__":
    main()
