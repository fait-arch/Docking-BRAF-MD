> En este repositorio solo están los scripts, `entrada/` y `sistema/`; `md_prueba/`, `analisis/` y
> `resultados_prueba/` se generan al correr los scripts.

# Pallidine (PubChem 24543347) + BRAF V600E: dinámica molecular

Mismo protocolo, scripts y receptor que `../Macaridine_94717379` (réplica de la sección 2.1.3 de
**RSC Adv. 2022, 12, 30181**): AmberTools para preparar el sistema y OpenMM como motor. Solo cambia el ligando.

## Entradas (`entrada/`, copia de `Simulacion/Data`)

- `Pallidine24543347_receptor.pdb`: **idéntico byte a byte** al receptor de Macaridine (BRAF V600E de
  SWISS-MODEL). Tras recorte 448–723 y PDB2PQR/PROPKA a pH 7.4, `sistema/receptor_amber.pdb` también es
  idéntico al de Macaridine (mismos HID/HIE, Glu600 como GLH).
- `Pallidine24543347_ligand.pdb`: SwissParam, 45 átomos, neutro. **No es una pose acoplada**: queda dentro de la
  proteína (54 pares a < 2.2 Å con Ser616, Asp576, Gly615...), a 16.8 Å del sitio ATP. Por eso se hace docking con Vina.

## Resultado de la prueba rápida (20 ps + 20 ps + 100 ps)

| | Pallidine | Macaridine |
|---|---|---|
| Vina, mejor pose | −8.07 kcal/mol (63 poses) | −7.38 kcal/mol |
| RMSD Cα, producción | 1.39 ± 0.08 Å | 1.37 ± 0.10 Å |
| RMSD ligando, producción | **3.89 ± 0.36 Å** (aún sube) | 1.15 ± 0.25 Å |
| Ligando enterrado | 88.5 % | 94.6 % |
| Puente de H principal | ninguno (< 5 %) | Thr529, 99 % |
| MM/GBSA ΔG (sin −TΔS) | **−34.4 ± 2.5 kcal/mol** | −27.4 ± 2.9 kcal/mol |
| ΔG por átomo pesado | −1.43 | −1.71 |
| Residuos dominantes | Phe583, Leu514, Val471 | Thr529, Lys483, Leu514 |

- Pallidine se queda en el sitio de ATP pero se mueve como bloque rígido (RMSD interno 0.54 Å): se desplaza
  ~3 Å y gira ~45° hacia la entrada del bolsillo (Asn580, Asn581, Lys591, Gln530, Ser536). Pierde el contacto
  con Thr529 (gatekeeper) y el contacto C=O···HN de la bisagra (Cys532) que tenía la pose de docking.
- MM/GBSA más favorable que Macaridine solo por tamaño (van der Waals −42.6 frente a −30.9 kcal/mol);
  electrostática −6.0 frente a −15.0.
- La amina terciaria queda a ~4.6 Å del carboxilato de Asp594. Se simuló **neutra**; a pH 7.4 probablemente
  está protonada, y la forma catiónica podría formar un puente salino. Es lo primero que conviene probar:
  `python 01_preparar_sistema.py ... --carga-ligando 1` (tras añadir el H al N en el PDB de entrada).
- Sin NaN. La ejecución se pausó durante el segmento 1 de producción; ese segmento se repitió completo
  desde `md_prueba/03_eq.xml`.
- Velocidad: 1.2 ns/día (NVT), 1.0 (equilibrado), 0.8 ns/día (producción) con 2 CPU.

## Cómo correrlo

```bash
micromamba create -f environment.yml && micromamba activate md-braf
cd scripts
python 01_preparar_sistema.py --receptor ../entrada/Pallidine24543347_receptor.pdb \
       --ligando ../entrada/Pallidine24543347_ligand.pdb --out ../sistema
python 02_dinamica.py --sistema ../sistema --out ../md_prueba \
       --t-calentamiento 20 --t-equilibrado 20 --t-produccion 100 --segmentos 4 --guardar-cada 1
python 03_rmsd.py --sistema ../sistema --md ../md_prueba --out ../resultados_prueba
python 04_analisis_reporte.py --base .. --md md_prueba --fig ../resultados_prueba --out ../analisis
# MM/GBSA (topologías y entrada en analisis/mmgbsa/)
cd ../analisis/mmgbsa && MMPBSA.py -O -i mmgbsa.in -o FINAL_MMGBSA.dat -do FINAL_DECOMP.dat \
   -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y ../../md_prueba/04_prod_00*.dcd \
   -eo energias_por_frame.csv -deo decomp_por_frame.csv
cd ../../scripts && python 05_mmgbsa_figura.py --mmgbsa ../analisis/mmgbsa --fig ../resultados_prueba
pymol -cq figuras_pymol.pml     # figuras 3D en resultados_prueba/
```

Nota: con la versión de `environment.yml` (python 3.11) conda puede elegir AmberTools 24.8, que construye una
caja octaédrica más pequeña que AmberTools 26 (a = 84.9 Å frente a 89.5 Å). Para reproducir este sistema, fijar
`ambertools=26.0`. PyMOL es incompatible con AmberTools 26 en el mismo entorno; aquí se usó un entorno aparte.

## Cambios en los scripts respecto a Macaridine

- `04_analisis_reporte.py`: el SMILES de la figura 2D se lee de `sistema/sistema_info.json` (antes estaba fijo
  el de Macaridine) y se dibuja con CoordGen (mejor para anillos puenteados).
- `03_rmsd.py`, `01_preparar_sistema.py`, `correr_50ns.sh`: solo el nombre del compuesto.
- `figuras_pymol.pml`: residuos del bolsillo de Pallidine y figura de dos paneles (pose inicial / final).

## Contenido

- `sistema/`: topología y coordenadas AMBER, parámetros del ligando (`LIG.mol2`, `LIG.frcmod`), poses de Vina.
- `md_prueba/`: minimización, calentamiento, equilibrado y 4 segmentos de producción (`.dcd`, `.log`, `.xml`, `.pdb`).
- `analisis/`: CSV de termodinámica, series, RMSF, contactos, puentes de H, poses de Vina, `resumen.json` y `mmgbsa/`.
- `resultados_prueba/`: `rmsd.png/csv`, figuras 3D (`fig_complejo_general.png`, `fig_bolsillo*.png`).

## Reporte

Sección 7 de `Docking/simulacion.tex` (y `simulacion.pdf`); figuras en `Docking/Img/Pallidine_24543347/`.
