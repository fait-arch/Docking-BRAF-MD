# Docking-BRAF-MD

*Docking* (AutoDock Vina) y dinámica molecular (AMBER ff19SB/GAFF + OpenMM) de compuestos naturales en el
sitio de ATP de la quinasa **BRAF V600E**, siguiendo el protocolo de la sección 2.1.3 de
*RSC Adv.* 2022, **12**, 30181 (Singh *et al.*).

Este repositorio contiene **solo el programa y los datos de entrada**. Las trayectorias, los análisis y las
figuras **no** se versionan: se generan al correr los scripts (ver `.gitignore`).

| Compuesto | PubChem CID | Carpeta |
|---|---|---|
| Macaridine | 94717379 | `simulaciones/Macaridine_94717379/` |
| Pallidine | 24543347 | `simulaciones/Pallidine_24543347/` |

## Estructura

```
Data/                              entradas originales (receptor SWISS-MODEL + ligando SwissParam)
environment.yml                    entorno conda/micromamba (md-braf)
simulaciones/<Compuesto_CID>/
    scripts/                       01 preparar · 02 dinámica · 03 RMSD · 04 análisis · 05 MM/GBSA
                                   benchmark.py · correr_50ns.sh · figuras_pymol.pml
    sistema/                       sistema AMBER ya preparado (topología, coordenadas, ligando GAFF, poses de Vina)
    entrada/                       (Pallidine) copia de las entradas de Data/
    README.md                      detalles y resultados de la prueba rápida de ese compuesto
```

`sistema/` es la **entrada de la dinámica**: lo genera `01_preparar_sistema.py` (recorte 448–723,
protonación a pH 7.4, *docking* con Vina, parámetros GAFF/AM1-BCC, caja TIP3P y NaCl 0.15 M). Se incluye para
que cualquier PC simule exactamente el mismo sistema, sin repetir el *docking*. Si se regenera, Vina puede
dar una pose algo distinta según la versión y el equipo.

## Instalación (una vez por PC)

```bash
# micromamba (si no está instalado)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)      # luego abrir una terminal nueva

git clone https://github.com/fait-arch/Docking-BRAF-MD.git
cd Docking-BRAF-MD
micromamba create -f environment.yml
micromamba activate md-braf
python -m openmm.testInstallation                   # con GPU NVIDIA debe aparecer CUDA sin error
```

Si CUDA da `CUDA_ERROR_UNSUPPORTED_PTX_VERSION`, el entorno trae un CUDA más nuevo que el driver. Fijarlo a
la versión que muestra `nvidia-smi`, por ejemplo `micromamba install -c conda-forge "cuda-version=13.2"`.

### Driver con CUDA 12.0–12.8 (por ejemplo, el HPC de CEDIA)

Si `nvidia-smi` muestra `CUDA Version` entre 12.0 y 12.8, OpenMM 8.6 no sirve: sus builds piden CUDA ≥ 12.9.
Usar el entorno alternativo, con OpenMM 8.2 (misma física; los scripts funcionan con ambas versiones):

```bash
micromamba create -f environment-cuda12.yml      # trae cuda-version=12.2; cambiarlo si nvidia-smi dice otro
micromamba activate md-braf
python -m openmm.testInstallation                # CUDA debe salir sin error
```

En el HPC de CEDIA (A100 en particiones MIG, driver 535 = CUDA 12.2), el escritorio interactivo
*NVIDIA DeepOps Desktop* funciona con: 1 GPU de 10 GB, 4 núcleos, 16 GB de RAM, 16 h e *idle timeout*
de 57600 s. Dentro de la sesión, lanzar la corrida con `nohup` para poder cerrar el visor:

```bash
cd simulaciones/<Compuesto_CID>/scripts
python 02_dinamica.py --sistema ../sistema --out ../md_prueba --etapas min
REUSAR_MIN=1 nohup bash correr_50ns.sh > corrida_50ns.log 2>&1 &
tail -f corrida_50ns.log
```

## Correr una simulación

Desde la carpeta `scripts/` del compuesto, con el entorno activo:

```bash
cd simulaciones/Pallidine_24543347/scripts

# velocidad esperada: minimizar una vez (unos minutos) y medir
python 02_dinamica.py --sistema ../sistema --out ../md_prueba --etapas min
python benchmark.py --sistema ../sistema --estado ../md_prueba/01_min.xml --pasos 2000

# 50 ns: minimización, 0.5 ns de calentamiento, 0.5 ns de equilibrado, 50 ns de producción,
# RMSD, análisis completo y MM/GBSA (200 fotogramas). Reanudable: si se corta, relanzar igual.
# REUSAR_MIN=1 aprovecha la minimización de md_prueba/ si ya existe.
REUSAR_MIN=1 systemd-inhibit --what=sleep:idle bash correr_50ns.sh

# otras duraciones:  T_PROD_PS=100000 NOMBRE=100ns bash correr_50ns.sh
```

Prueba rápida (20 ps + 20 ps + 100 ps, la del reporte):

```bash
python 02_dinamica.py --sistema ../sistema --out ../md_prueba \
       --t-calentamiento 20 --t-equilibrado 20 --t-produccion 100 --segmentos 4 --guardar-cada 1
python 03_rmsd.py --sistema ../sistema --md ../md_prueba --out ../resultados_prueba
```

Los resultados quedan en `md_*/`, `analisis_*/` y `resultados_*/` dentro de la carpeta del compuesto,
ignorados por git. Velocidad de referencia con ~49 000 átomos: ~0.8 ns/día con 2 CPU y ~135 ns/día con una
RTX 3060 de portátil (CUDA). Con GPU, 50 ns tardan unas 9 h por compuesto, más el análisis.

## Añadir un compuesto

1. Poner `<Nombre><CID>_ligand.pdb` y `<Nombre><CID>_receptor.pdb` en `Data/`.
2. Copiar `simulaciones/Pallidine_24543347/` como `simulaciones/<Nombre>_<CID>/` y borrar su `sistema/` y
   `entrada/`. Sus scripts leen el SMILES del ligando de `sistema/sistema_info.json`; los de Macaridine lo tienen fijo.
3. Cambiar el nombre del compuesto en `03_rmsd.py` (`--titulo`) y en `correr_50ns.sh`.
4. Preparar el sistema:
   ```bash
   python 01_preparar_sistema.py --receptor ../../../Data/<Nombre><CID>_receptor.pdb \
          --ligando ../../../Data/<Nombre><CID>_ligand.pdb --out ../sistema
   ```
   Revisar la carga del ligando (`--carga-ligando`): por ejemplo, una amina protonada a pH 7.4 es +1.

## Trabajar desde otra PC

```bash
git clone https://github.com/fait-arch/Docking-BRAF-MD.git    # o: gh repo clone fait-arch/Docking-BRAF-MD
# ... cambios ...
git add -A && git commit -m "describe el cambio" && git push
```

Para hacer `push` hay que autenticarse una vez en cada PC, con `gh auth login` o con una clave SSH. Antes de
empezar a trabajar en una PC, `git pull` trae los cambios hechos en la otra.

## Notas

- Las versiones usadas en el reporte son AmberTools 26, OpenMM 8.6.1, AutoDock Vina 1.2.7, RDKit 2026.03,
  ParmEd 4.3.1, PDB2PQR 3.6–3.8 y MDAnalysis 2.10–2.11.
- `correr_50ns.sh` genera las topologías de MM/GBSA (radios mbondi2) si no encuentra `analisis/mmgbsa/`.
- El reporte LaTeX (`simulacion.tex`) y sus figuras no forman parte de este repositorio.
