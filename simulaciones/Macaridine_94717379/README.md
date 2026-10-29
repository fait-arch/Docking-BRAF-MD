> En este repositorio solo están los scripts y `sistema/`; `md_prueba/`, `analisis/` y `resultados_prueba/` se
> generan al correr los scripts. Los datos de entrada están en `../../Data/`.

# Macaridine (PubChem 94717379) + BRAF V600E: dinámica molecular

Réplica del protocolo de MD de **RSC Adv. 2022, 12, 30181 (d2ra05751d), sección 2.1.3**,
con AmberTools para preparar el sistema y OpenMM como motor. El artículo usa Amber,
así que se mantienen el mismo campo de fuerza y los mismos parámetros. El módulo
`automd` de dockingML depende de GROMACS y su generador de topologías no tiene
ff19SB, la caja octaédrica ni el cálculo de iones SLTCAP.

## Parámetros: artículo vs. prueba rápida

| | Artículo | Esta carpeta |
|---|---|---|
| Proteína | AMBER ff19SB, dominio quinasa BRAF V600E (PDB 4XV2) | ff19SB, residuos 448–723 del modelo SWISS-MODEL (V600E) |
| Protonación | PDB2PQR, pH fisiológico | PDB2PQR + PROPKA, pH 7.4 |
| Ligando | GAFF (Antechamber, AMBER 20) | GAFF + AM1-BCC (Antechamber, AmberTools 26) |
| Pose inicial | Docking con AutoDock Vina, caja 20 Å, exhaustiveness = 100, 100 modos | Igual (caja centrada en el bolsillo ATP) |
| Agua / caja | TIP3P, octaedro truncado, 12 Å | Igual |
| Iones | Neutralización + NaCl 0.15 M (Machado & Pantano, SLTCAP) | Igual: carga +8 → 36 Na⁺ / 44 Cl⁻ |
| Minimización | 10 000 ciclos (4000 SD + 6000 CG), restricción 10 kcal/mol/Å² | 10 000 iteraciones L-BFGS (OpenMM), misma restricción sobre átomos pesados |
| Calentamiento | 0 → 310 K en 500 ps, sin restricciones | 0 → 310 K en 31 escalones · **prueba: 20 ps** |
| Equilibrado | 500 ps | NPT 310 K / 1 atm · **prueba: 20 ps** |
| Producción | 1 µs NPT (300 ns para dabrafenib) | NPT · **prueba: 100 ps** |
| Termostato / barostato | Langevin, γ = 1 ps⁻¹ · 1 atm | Langevin (LangevinMiddle), γ = 1 ps⁻¹ · Monte Carlo 1 atm |
| Paso / restricciones | 2 fs, SHAKE en enlaces con H | 2 fs, enlaces con H restringidos |
| No enlazadas | Corte 11 Å, PME, PBC | Igual |

## Avisos sobre los datos de entrada (`Simulacion/Data`)

1. **`ligand.pdb` no está acoplado a la proteína.** Son las coordenadas de SwissParam
   (centradas en el origen) y caen *dentro* de la proteína: 71 contactos a < 2.2 Å con el
   esqueleto de Lys578, Gly615 y Ser616, a ~17 Å del sitio ATP. Por eso se repite el docking
   con Vina, como en el artículo. La simulación del Colab (`SimulacionesRespados/94717379_*.ipynb`)
   partió de esta pose solapada.
2. **`receptor.pdb` es BRAF completa (766 residuos)**, con dominios desordenados. Con 12 Å de
   agua el sistema pasa de 400 000 átomos. El artículo simula solo el dominio quinasa, así que
   aquí se recorta a 448–723 (`--residuos`). Para usar la proteína completa: `--residuos 1-766`.
3. PROPKA asigna **Glu600 protonado (GLH)** en este modelo. Las histidinas quedan como
   HID/HIE según `sistema/sistema_info.json`.

## Resultado de la prueba rápida (`resultados_prueba/rmsd.png`)

Referencia: el complejo minimizado. Producción de 100 ps a 310.5 ± 1.5 K y 1.009 g/mL.

| | Calentamiento | Equilibrado | Producción |
|---|---|---|---|
| RMSD Cα proteína (Å) | 0.96 ± 0.06 | 1.10 ± 0.04 | **1.37 ± 0.10** (máx. 1.57) |
| RMSD ligando, átomos pesados (Å) | 1.26 ± 0.33 | 1.86 ± 0.34 | **1.15 ± 0.25** (máx. 1.79) |

- El ligando se mantiene en el bolsillo ATP todo el tiempo. Contactos al final: Val471, Ala481,
  Lys483, Leu514, Thr529 (gatekeeper), Trp531 y Cys532 (bisagra), Asn581, Phe583, Gly593 y Asp594 (DFG).
- El RMSD de la proteína todavía sube (1.2 → 1.55 Å): con 100 ps no hay convergencia, como era de
  esperar. Lo que más se mueve es el lazo de activación (608–614), la región 658–660 y los extremos
  recortados (448, 722–723).
- Valores comprobados con un segundo programa (mdtraj): coinciden con MDAnalysis.
- Velocidad aquí: ~0.8 ns/día en NPT con 2 CPU. Con GPU (OpenMM CUDA) debería ser de cientos de
  ns/día; `benchmark.py` lo mide.

**Incidencia:** el segundo segmento de producción dio NaN a los ~28 ps. El estado previo estaba limpio
(ningún par de átomos a < 1.3 Å, temperatura y densidad estables). Se repitió desde
`04_prod_001.xml` y llegó a 100 ps sin más problemas. `02_dinamica.py` ahora reintenta un segmento
automáticamente si aparece un NaN y lo registra en `md_prueba/incidencias.log`. Si se repite en la
corrida larga, vale la pena revisarlo con más detalle.

## Cómo correrlo

```bash
# entorno (una vez)
micromamba create -f environment.yml      # o: conda env create -f environment.yml
micromamba activate md-braf

cd scripts
# 1) sistema: recorte, protonación, docking Vina, GAFF, tleap
python 01_preparar_sistema.py --receptor ../../../Data/Macaridine94717379_receptor.pdb \
       --ligando ../../../Data/Macaridine94717379_ligand.pdb --out ../sistema

# 2a) prueba rápida (lo que hay en md_prueba/)
python 02_dinamica.py --sistema ../sistema --out ../md_prueba \
       --t-calentamiento 20 --t-equilibrado 20 --t-produccion 100 --segmentos 4 --guardar-cada 1

# 2b) protocolo completo del artículo (500 ps + 500 ps + 1 µs en 20 segmentos de 50 ns).
#     Usa CUDA automáticamente si hay GPU NVIDIA; es reanudable.
python 02_dinamica.py --sistema ../sistema --out ../md

# 3) RMSD
python 03_rmsd.py --sistema ../sistema --md ../md_prueba --out ../resultados_prueba

# velocidad esperada en tu equipo
python benchmark.py --sistema ../sistema --estado ../md_prueba/01_min.xml
```

Las trayectorias guardan solo proteína + ligando (se leen con `sistema/complejo_seco.prmtop`);
1 µs con un fotograma cada 10 ps ocupa unos 5.4 GB. Con `--dcd-completo` se guarda también el agua.
En esta prueba, las trayectorias de calentamiento y equilibrado (`02_cal.dcd`, `03_eq.dcd`) se
escribieron con el sistema completo. `03_rmsd.py` detecta cada caso solo.

## Contenido

- `scripts/`: preparación (01), dinámica (02), RMSD (03), análisis del reporte (04), MM/GBSA (05),
  `benchmark.py` y `figuras_pymol.pml`.
- `sistema/`: topología y coordenadas AMBER (`sistema.prmtop/.inpcrd/.pdb`), parámetros del ligando
  (`LIG.mol2`, `LIG.frcmod`), poses de Vina, entradas de tleap, `sistema_info.json`.
- `md_prueba/`: minimización, calentamiento, equilibrado y producción (`.dcd`, `.log`, `.xml`, `.pdb`).
- `resultados_prueba/`: `rmsd.png`, `rmsd.csv`, `rmsd_resumen.csv`.
- `analisis/`: CSV de termodinámica, series estructurales, RMSF, contactos, puentes de H, poses de Vina,
  `resumen.json` y `mmgbsa/` (MM/GBSA GB-OBC2 con descomposición por residuo).

## Simulación de 50 ns

`scripts/correr_50ns.sh` hace todo de una vez con los parámetros del artículo: minimización,
calentamiento de 500 ps, equilibrado de 500 ps, **producción de 50 ns** (10 segmentos de 5 ns,
reanudables), RMSD, análisis completo y MM/GBSA (200 fotogramas). Necesita GPU: en CPU, 50 ns
tardarían semanas (aquí se midió ~0.8 ns/día).

```bash
micromamba create -f environment.yml && micromamba activate md-braf
python -m openmm.testInstallation          # debe aparecer CUDA
cd scripts
bash correr_50ns.sh                         # resultados en ../md_50ns, ../analisis_50ns, ../resultados_50ns
```

Otras duraciones: `T_PROD_PS=100000 NOMBRE=100ns bash correr_50ns.sh`. Para reutilizar la minimización
de la prueba (el sistema es el mismo): `REUSAR_MIN=1`. Si se interrumpe, se relanza igual y continúa
desde el último segmento terminado.

## Reporte

El reporte completo, con teoría, parámetros, gráficas y conclusión, está en
`Docking/simulacion.tex` (y `simulacion.pdf`); sus figuras están en `Docking/Img/Macaridine_94717379/`.
Se regeneran con:

```bash
cd scripts
python 04_analisis_reporte.py --base .. --fig ../resultados_prueba --out ../analisis
# MM/GBSA (AmberTools): ver analisis/mmgbsa/mmgbsa.in; luego
python 05_mmgbsa_figura.py --mmgbsa ../analisis/mmgbsa --fig ../resultados_prueba
pymol -cq figuras_pymol.pml        # figuras 3D (quedan en resultados_prueba/)
```

## Versión con machine learning (ML/MM)

`correr_ml.sh` repite el mismo protocolo con la energía interna del ligando calculada por la red neuronal
ANI-2x (`potencial_ml.py`), corre 20 ns (ampliables a 50 ns con `T_PROD_PS=50000`) y compara con la corrida
clásica de `md_50ns/` (`06_comparar_clasico_ml.py`). Instalación, tiempos e interpretación en `../../ML.md`.

```bash
nohup bash correr_ml.sh > corrida_ml_20ns.log 2>&1 &
T_PROD_PS=50000 nohup bash correr_ml.sh > corrida_ml_50ns.log 2>&1 &
```
