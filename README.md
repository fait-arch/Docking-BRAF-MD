# Docking-BRAF-MD

*Docking* (AutoDock Vina) y dinámica molecular (AMBER ff19SB/GAFF + OpenMM) de compuestos naturales en el
sitio de ATP de la quinasa **BRAF V600E**, siguiendo el protocolo de la sección 2.1.3 de
*RSC Adv.* 2022, **12**, 30181 (Singh *et al.*).

El repositorio contiene **el programa, los datos de entrada y el reporte**. Las trayectorias y los análisis
pesados no se versionan: se generan al correr los pipelines (ver `.gitignore`).

| Sistema | Clave para la CLI | PubChem CID | Carga | 10 ns | Nota |
|---|---|---|---|---|---|
| Dabrafenib | `dabrafenib` | 44462760 | 0 | ✅ | **control positivo** |
| Macaridine | `macaridine` | 94717379 | 0 | ✅ | |
| Notholaenic acid | `notholaenic_acid` | 78213884 | 0 | ✅ | forma neutra |
| Notholaenic acid anión | `notholaenic_acid_anion` | 78213884 | −1 | — | sin pKa medido |
| Pallidine | `pallidine` | 24543347 | 0 | ✅ | forma neutra |
| Pallidine protonada | `pallidine_protonada` | 24543347 | +1 | — | |
| Quinine | `quinine` | 72243208 | 0 | ✅ | protonación incorrecta |
| Quinine protonada | `quinine_protonada` | 72243208 | +1 | — | **la válida a pH 7.4** |

### Estados de protonación — leer antes de interpretar nada

A pH 7.4 varios de estos compuestos **no están neutros**, y la carga no es un número que se pueda cambiar
en el YAML: `mol_ligando()` deduce los enlaces de la geometría 3D contando los hidrógenos del archivo, así
que el PDB tiene que llevar el protón puesto o quitado. Cambiarla obliga además a rehacer las cargas
AM1-BCC, los contraiones de la caja y el *docking*, porque Vina puntúa la electrostática. Es decir: hay que
volver a correr la etapa `sistema` completa.

| Compuesto | Estado a pH 7.4 | Base |
|---|---|---|
| Quinine | **+1**, 91 % | pKa 8.4 — Roy L *et al.*, *Br J Clin Pharmacol* 2002;54:604-609 (PMID 12492607) |
| Notholaenic acid | **−1** (inferido) | carboxilo confirmado — Cioffi G *et al.*, *Molecules* 2011;16:2527-2541 (PMID 21415834). **Sin pKa experimental publicado**: se corren las dos formas |
| Pallidine | **+1** | amina terciaria |
| Dabrafenib, Macaridine | 0 | |

Los PDB cargados los genera `generar_ligandos.py` desde el SMILES, y cada uno se verifica releyéndolo con
el mismo `mol_ligando()` del pipeline y comparando su InChIKey. Las formas neutras se conservan como
entradas aparte: son las que tienen los 10 ns corridos y hacen falta para la comparación.

---

## Cómo se ejecuta: una sola CLI, un parámetro por PC

Todo se lanza y se consulta con `dmlrun`. **Lo único que cambia entre máquinas es `--servidor`**: el gestor
de colas, los recursos, el entorno y la plataforma de OpenMM salen de `config/servidores.yaml`, así que no
hay que editar código ni scripts al cambiar de equipo.

```bash
dmlrun servidores      # perfiles de PC disponibles
dmlrun compuestos      # sistemas registrados y qué corridas tiene cada uno
dmlrun config --servidor hpc    # comprueba, en esta máquina, que exista lo que el perfil pide
dmlrun run             # lanza
dmlrun status          # estado y progreso real
dmlrun logs --id <id>  # log de una corrida
dmlrun reporte         # informe HTML, o secciones LaTeX dentro del paper
```

### Perfiles de PC

| `--servidor` | Dónde | Cola | Recursos |
|---|---|---|---|
| `hpc` | HPC CEDIA | Slurm, partición `interactive` | 1 GPU MIG 10 GB, 4 núcleos, 16 GB, 16 h |
| `hpc-cpu` | HPC CEDIA | Slurm | 16 núcleos, 48 GB, 24 h, sin GPU |
| `linux-gpu` | Servidor Linux externo | ninguna (`setsid`) | GPU NVIDIA |
| `linux-cpu` | Servidor Linux externo | ninguna (`setsid`) | solo CPU |
| `portatil` | Portátil con RTX 3060 (micromamba) | ninguna | 16-30 ns/día medidos |
| `actual` | Aquí mismo, entorno ya activo | ninguna | pruebas cortas y figuras |

`--servidor auto` elige solo: respeta `$DMLRUN_SERVIDOR`, luego un perfil que se llame como el *hostname*,
luego el primero de Slurm si existe `sbatch`, y si no uno que ejecute en local.

**Añadir una PC** es copiar un bloque en `config/servidores.yaml` y ya se puede usar su nombre:

```yaml
  mi-servidor:
    descripcion: "Servidor del laboratorio"
    scheduler: local            # local | slurm
    cpus: 16
    memoria_gb: 64
    gpus: 1
    plataforma_openmm: CUDA
    conda_env: md-braf
    gestor_entorno: micromamba
    pre_comandos:
      - 'eval "$(micromamba shell hook -s bash)"'
```

---

## Instalación (una vez por PC)

```bash
# micromamba (si no está instalado)
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)      # luego abrir una terminal nueva

git clone https://github.com/fait-arch/Docking-BRAF-MD.git
cd Docking-BRAF-MD
micromamba create -f environment.yml
micromamba activate md-braf
pip install -e .                                     # deja disponible el comando dmlrun
python -m openmm.testInstallation                    # con GPU NVIDIA debe aparecer CUDA sin error

dmlrun config --servidor <perfil-de-esta-pc>         # comprueba el entorno antes de lanzar nada
```

Sin `pip install -e .` todo funciona igual escribiendo `python -m dmlrun.cli` en lugar de `dmlrun`.

Si CUDA da `CUDA_ERROR_UNSUPPORTED_PTX_VERSION`, el entorno trae un CUDA más nuevo que el driver. Fijarlo a
la versión que muestra `nvidia-smi`, por ejemplo `micromamba install -c conda-forge "cuda-version=13.2"`.

### Driver con CUDA 12.0–12.8 (por ejemplo, el HPC de CEDIA)

Si `nvidia-smi` muestra `CUDA Version` entre 12.0 y 12.8, OpenMM 8.6 no sirve: sus builds piden CUDA ≥ 12.9.
Usar el entorno alternativo, con OpenMM 8.2 (misma física; los scripts funcionan con ambas versiones):

```bash
micromamba create -f environment-cuda12.yml      # trae cuda-version=12.2; cambiarlo si nvidia-smi dice otro
micromamba activate md-braf
pip install -e .
python -m openmm.testInstallation                # CUDA debe salir sin error
```

En el HPC de CEDIA (A100 en particiones MIG, driver 535 = CUDA 12.2) el escritorio interactivo
*NVIDIA DeepOps Desktop* funciona con 1 GPU de 10 GB, 4 núcleos, 16 GB de RAM, 16 h e *idle timeout* de
57600 s. Es exactamente lo que pide el perfil `hpc`.

---

## Las entradas

Hay tres niveles de entrada, de lo más general a lo más específico:

| Qué | Dónde | Quién lo usa |
|---|---|---|
| PDB originales del receptor y del ligando | `Data/<Nombre><CID>_{receptor,ligand}.pdb` | la etapa `sistema` |
| Metadatos del compuesto (SMILES, carga, título, bolsillo) | `config/compuestos.yaml` | toda la CLI |
| Sistema AMBER ya preparado | `simulaciones/<Sistema>/sistema/` | la dinámica |

`sistema/` es **la entrada real de la dinámica**: lo genera `01_preparar_sistema.py` (recorte 448–723,
protonación a pH 7.4, *docking* con Vina, parámetros GAFF/AM1-BCC, caja TIP3P y NaCl 0.15 M) y se versiona
para que cualquier PC simule exactamente el mismo sistema sin repetir el *docking*. Si se regenera, Vina
puede dar una pose algo distinta según la versión y el equipo.

Lo que distingue a un sistema de otro está en `config/compuestos.yaml`, no en los scripts: el SMILES, la
**carga del ligando** (0 para la forma neutra, +1 para la protonada), el título de las figuras y los residuos
que definen el centro del bolsillo. Los scripts de `simulaciones/scripts/` son los mismos para todos.

> **Pendiente.** Los tres sistemas con carga (`quinine_protonada`, `notholaenic_acid_anion`,
> `pallidine_protonada`) todavía no tienen su `sistema/` preparado. Sus PDB ya están en `Data/` y
> verificados; falta correr la etapa:
> ```bash
> ENVIAR=1 ./correr_hpc_cedia.sh prueba quinine_protonada,notholaenic_acid_anion,pallidine_protonada
> ```
> Hasta noviembre de 2026 `pallidine_protonada` era inejecutable: declaraba carga +1 apuntando al PDB
> **neutro** de pallidine, y `DetermineBondOrders` no puede colocar una carga sobre una molécula con los
> hidrógenos de la forma neutra. Ya tiene su propio PDB.

### Añadir un compuesto

1. Poner `<Nombre><CID>_ligand.pdb` y `<Nombre><CID>_receptor.pdb` en `Data/`.
2. Añadir un bloque en `config/compuestos.yaml`:
   ```yaml
     mi_compuesto:
       nombre: MiCompuesto
       pubchem: 12345678
       directorio: simulaciones/MiCompuesto_12345678
       titulo: "MiCompuesto (12345678) – BRAF V600E"
       receptor: MiCompuesto12345678_receptor.pdb
       ligando: MiCompuesto12345678_ligand.pdb
       smiles: "..."
       carga_ligando: 0          # +1 si tiene una amina protonada a pH 7.4
   ```
3. Preparar el sistema y lanzar:
   ```bash
   dmlrun run --pipeline resultados --compuesto mi_compuesto --tiempo 50 \
          --etapas sistema --servidor actual
   dmlrun run --pipeline resultados --compuesto mi_compuesto --tiempo 50 --servidor hpc
   ```

No hay que copiar scripts ni editar ninguno: antes había que cambiar el nombre a mano en `03_rmsd.py` y en
`correr_50ns.sh`, y ya no.

---

## Los dos pipelines

### `--pipeline resultados` — dinámica tradicional

Etapas: `sistema`, `min`, `cal`, `eq`, `prod`, `rmsd`, `analisis`, `mmgbsa`.

```bash
# todos los compuestos registrados, 50 ns, en la HPC
dmlrun run --pipeline resultados --compuesto todos --tiempo 50 --servidor hpc

# ver el plan y el sbatch sin enviar nada
dmlrun run --pipeline resultados --compuesto todos --tiempo 50 --servidor hpc --dry-run
dmlrun run --pipeline resultados --compuesto macaridine --tiempo 50 --servidor hpc \
       --dry-run --mostrar-script

# trabajos más cortos: como máximo 2 segmentos de producción por envío
dmlrun run --pipeline resultados --compuesto macaridine --tiempo 50 \
       --etapas prod --max-segmentos 2 --servidor hpc

# análisis y MM/GBSA sin gastar cuota de GPU
dmlrun run --pipeline resultados --compuesto todos --tiempo 50 \
       --etapas rmsd,analisis,mmgbsa --servidor hpc-cpu

# prueba rápida antes de ocupar la cola (100 ps, mismo protocolo)
dmlrun run --pipeline resultados --compuesto macaridine --tiempo prueba --servidor actual
```

El tiempo de simulación es un parámetro, no una carpeta fija: `--tiempo 10,20,50` lanza las tres y cada una
queda en `runs/<tiempo>/` sin pisar a las otras.

### `--pipeline ml` — docking y rescoring con aprendizaje automático

Etapas: `centro`, `docking`, `complejos`, `features`, `rescore`.

Necesita además el paquete `dockingML`, que aporta los comandos `dml-dock`, `dml-features`, `dml-train` y
`dml-rescore`:

```bash
pip install -e /ruta/a/dockingML
```

```bash
# sobre las poses acopladas: no depende de la dinámica
dmlrun run --pipeline ml --compuesto todos --tiempo 50 --fuente docking --servidor hpc-cpu

# sobre fotogramas de la trayectoria de 50 ns
dmlrun run --pipeline ml --compuesto todos --tiempo 50 --fuente md --frames-ml 20 --servidor hpc-cpu

# con un modelo entrenado
dmlrun run --pipeline ml --compuesto todos --tiempo 50 \
       --fuente docking --modelo modelos/rescorer.pkl --servidor hpc-cpu
```

El rescoring necesita un modelo entrenado con afinidades conocidas frente a BRAF (por ejemplo de ChEMBL);
sin él, el pipeline deja los descriptores calculados en `features.csv` y no puntúa.

---

## Dónde se guardan los resultados

Cada ejecución estrena su propia carpeta, fechada, dentro de `resultados/`. **Nada se sobrescribe y nada
va a GitHub** (`/resultados/` está en `.gitignore`):

```
resultados/
├── 20261105-1702/                  una carpeta por TANDA de lanzamiento
│   └── 10nsMacaridine/
├── 20261106-0107/                  los tres que se lanzaron juntos
│   ├── 10nsDabrafenib/
│   │   ├── md/          trayectorias .dcd y estados .xml (el 99 % del peso)
│   │   ├── resultados/  las 10 figuras, en PNG y PDF
│   │   ├── analisis/    resumen.json, CSV y mmgbsa/FINAL_MMGBSA.dat
│   │   └── estado/      log, script enviado, manifest.json
│   ├── 10nsQuinine/
│   └── 10nsNotholaenicAcid/
└── 20261108-0439/                  la tanda siguiente, sin mezclarse
    ├── 0.1nsQuinineProtonada/
    ├── 0.1nsNotholaenicAcidAnion/
    └── 0.1nsPallidineProtonada/
```

El sello `AAAAMMDD-HHMM` se calcula **una vez por invocación de la CLI**, no por
compuesto: así una tanda de cinco cae entera en la misma carpeta aunque el envío cruce un minuto.

Las corridas de disposiciones anteriores (planas en la raíz, o con la fecha pegada al nombre) se siguen
reconociendo, así que nada queda huérfano. Para agruparlas:

```bash
python3 migrar_resultados.py            # enseña el plan, no mueve nada
python3 migrar_resultados.py --aplicar
```

Toma la fecha del `run_id` de cada `manifest.json`, que es el momento real del envío.

Cuándo reanuda y cuándo estrena carpeta:

| Situación | Qué hace |
|---|---|
| Relanzas un compuesto ya corrido | **reanuda** — si Slurm cortó por tiempo, no se pierden horas |
| Compuesto o tiempo nuevos | carpeta nueva con fecha `AAAAMMDD-HHMM` |
| `--nueva-corrida` | carpeta nueva aunque haya una previa |
| **Rehiciste la topología** | **carpeta nueva automática** |

Lo último es una salvaguarda, no un adorno: compara la fecha de `sistema.prmtop` con la del `01_min.xml`
de la corrida. Sin ella, reanudar tras reconstruir un sistema mezclaría una trayectoria vieja con una
topología nueva, el pipeline daría por buenos los `.xml` que encuentra y el resultado sería un híbrido
silencioso.

Para empaquetar y comparar, sin trayectorias:

```bash
python3 exportar_resultados.py 10     # exportacion/resumen_10ns.csv + 10ns.tar.gz
python3 figura_rmsd.py 10             # exportacion/figuras/ (PNG 200 dpi y PDF vectorial)
```

Ninguno de los dos necesita el entorno activado: solo usan la librería estándar (el de figuras, además,
matplotlib).

---

## Seguimiento

```bash
dmlrun status                       # últimas corridas, con el progreso real
dmlrun status --compuesto macaridine
dmlrun logs --id <id> -f            # log en vivo
dmlrun cancel --id <id>
dmlrun comparar --compuesto todos   # qué hay calculado de cada sistema y tiempo
```

La columna **progreso** se lee del disco (`producción 7/10` cuenta los segmentos terminados), no del gestor
de colas, así que sigue siendo correcta aunque Slurm haya matado el trabajo por tiempo.

Los envíos son **reanudables**: `dmlrun run` omite las etapas que ya tienen salida. Si una corrida se corta,
se relanza el mismo comando y continúa donde estaba. `--forzar` repite todo.

---

## El reporte

El reporte LaTeX vive ahora en `paper/simulacion.tex`, con sus figuras en `paper/Img/`.

```bash
cd paper && latexmk -pdf simulacion.tex
```

Las secciones con cifras se generan de los datos, no se copian a mano:

```bash
dmlrun reporte --formato tex          # escribe dentro de paper/simulacion.tex
dmlrun reporte                        # informe HTML autocontenido en reportes/
```

`--formato tex` **solo reemplaza lo que hay entre sus propias marcas**
`% <<<dmlrun:...>>> ... % <<<dmlrun:fin:...>>>`, y guarda una copia `.bak` antes de escribir: la redacción
del documento no se toca. Si las marcas no existen, se insertan una vez antes del apéndice; moverlas de
sitio cambia dónde aparecen las secciones. Genera dos bloques:

- **`datos_50ns`** — tablas de estabilidad, energías y descomposición del MM/GBSA, leídas de
  `resumen.json` y `FINAL_MMGBSA.dat`, de modo que no pueden desfasarse del análisis.
- **`comparacion_md_ml`** — la comparación entre la vía tradicional (docking + 50 ns + MM/GBSA) y la de
  aprendizaje automático (docking + rescoring) sobre los mismos compuestos.

El informe HTML lleva las imágenes incrustadas y cada figura con qué muestra y qué hay que mirar en ella;
sirve para pasárselo a alguien que no lanzó la simulación.

---

## Reproducibilidad entre PC

Cada envío escribe `runs/<tiempo>/estado/manifest.json` con la **huella** de los parámetros científicos
(compuesto, SMILES, carga, pH, residuos, tiempo, segmentos y semilla), el commit de git y las versiones de
OpenMM, MDAnalysis, RDKit y demás. Misma huella en dos PC = resultados comparables; cambiar el tiempo de
simulación cambia la huella a propósito.

La semilla del integrador es fija (2022) y entra en la huella. Entre plataformas distintas (CUDA frente a
CPU, u OpenMM 8.2 frente a 8.6) las trayectorias no son idénticas bit a bit; el manifiesto deja constancia
de con qué se calculó cada una.

### Trabajar desde otra PC

```bash
git clone https://github.com/fait-arch/Docking-BRAF-MD.git    # o: gh repo clone fait-arch/Docking-BRAF-MD
git pull                                                      # antes de empezar, en cada PC
# ... cambios ...
git add -A && git commit -m "describe el cambio" && git push
```

Para hacer `push` hay que autenticarse una vez en cada PC, con `gh auth login` o con una clave SSH.

---

## Velocidad de referencia

Medido en noviembre de 2026 sobre este sistema (~49 000 átomos), no estimado:

| Máquina | Producción | 10 ns completos |
|---|---|---|
| HPC CEDIA, A100 MIG | **260 ns/día** | ~2 h con análisis y MM/GBSA |
| HPC CEDIA, trozo MIG pequeño | 125 ns/día | ~3 h 45 |
| Portátil RTX 3060 a 44 W, 2 trabajos a la vez | **15.8 ns/día** | ~17 h |
| CPU (2 núcleos) | ~0.8 ns/día | semanas |

> El portátil daba 15.8 ns/día con la GPU limitada a **44 W de los 115 disponibles** (`SW Power Cap` activo,
> relojes a 1522 de 2100 MHz, perfil del sistema en `balanced`). Se corrige en caliente con
> `powerprofilesctl set performance`. Y lanzar varios compuestos a la vez **no acelera nada**: con la GPU
> al 100 % se reparten el tiempo por rodajas.

Para medirlo en una PC nueva:

```bash
dmlrun run --pipeline resultados --compuesto macaridine --tiempo prueba \
       --etapas min --servidor <perfil>
python simulaciones/scripts/benchmark.py --sistema simulaciones/Macaridine_94717379/sistema \
       --estado <ruta al 01_min.xml> --pasos 2000
```

---

## Problemas frecuentes

**`sbatch: Requested node configuration is not available`.** Ningún nodo de la partición existe con esa
combinación; `dmlrun` lo explica al fallar. Revisar en este orden:

```bash
scontrol show node compute-0-0 | grep -i gres    # nombre real del GRES de GPU
sinfo -p interactive -o "%P %c %m %G"            # CPU, memoria y GPU por nodo
scontrol show partition interactive              # límite de tiempo
```

En nodos con GPU particionada (MIG) el GRES no es `gpu:1` sino algo como `gpu:2g.10gb`. Se ajusta en
`config/servidores.yaml`.

**`micromamba: command not found`, o `/micromamba: No such file or directory`.** El binario vive en
`~/.local/bin` y no siempre está en el `PATH`; además hace falta `MAMBA_ROOT_PREFIX`. En una línea:

```bash
export MAMBA_ROOT_PREFIX="$HOME/micromamba" PATH="$HOME/.local/bin:$PATH" && eval "$(micromamba shell hook -s bash)" && micromamba activate md-braf
```

`correr_hpc_cedia.sh` ya lo hace solo. En un trabajo de Slurm va en `pre_comandos` del perfil.

**Los trabajos se quedan en `PENDING` con `QOSMaxJobsPerUserLimit`.** La QOS `normal` de CEDIA permite **un
solo trabajo corriendo a la vez** (`sacctmgr -p show assoc user=$USER`), así que los compuestos van en serie.
Y una sesión de escritorio de Open OnDemand **es un trabajo de Slurm**: ocupa ese cupo mientras esté viva.
Para liberarlo sin perder la terminal, entra por **Clusters → Shell Access** (no consume cola) y cancela la
sesión de escritorio. Pedir una QOS con más cupo a soporte es la mejora de rendimiento más barata del
proyecto.

**No tengo SSH al clúster.** No hace falta: OnDemand trae explorador de archivos web en
`https://hpc.cedia.edu.ec/pun/sys/dashboard/files/fs/<ruta absoluta>`, con previsualización de PNG y botón
de descarga, y una terminal en **Clusters → Shell Access**.

**La terminal del escritorio VNC parte los comandos de varias líneas al pegarlos.** Usa comandos de una
sola línea, o guiones del repositorio.

**MM/GBSA.** El número de fotogramas y el intervalo se calculan contando los que realmente hay en los `.dcd`,
así que funciona igual si la producción quedó a medias.

---

## Estructura

```
Data/                              entradas originales (receptor SWISS-MODEL + ligando SwissParam)
config/
    servidores.yaml                perfiles de PC: cola, recursos, entorno
    compuestos.yaml                sistemas: SMILES, carga, bolsillo, rutas
dmlrun/                            la CLI de orquestación
simulaciones/
    scripts/                       01 preparar · 02 dinámica · 03 RMSD · 04 análisis
                                   05 MM/GBSA · 06 entrada de ML · benchmark
    <Sistema>/sistema/             sistema AMBER preparado (entrada de la dinámica)
resultados/<sesion>/<tiempo>ns<Nombre>/  corridas, por tanda (NO versionadas)
exportacion/                       CSV, paquetes y figuras combinadas (NO versionada)
migrar_resultados.py               agrupa corridas sueltas en carpetas por sesión
trayectorias_50ns/<Sistema>/       DCD de producción y topología seca
correr_hpc_cedia.sh                lanzador para CEDIA, con comprobaciones previas
generar_ligandos.py                PDB de los estados de protonación, desde SMILES
exportar_resultados.py             resumen comparativo en CSV + paquete sin trayectorias
figura_rmsd.py                     figura de RMSD combinada de todos los compuestos
paper/
    simulacion.tex                 el reporte
    Img/                           figuras del reporte
docs/EJECUCION.md                  guía larga de la CLI
```

---

## Estado del proyecto y advertencias (nov 2026)

Lo que hay que tener presente antes de interpretar o continuar nada.

**El ranking por MM/GBSA total todavía no es defendible.** Los 10 ns de los cinco compuestos terminaron bien,
pero el control positivo queda **tercero**:

| Compuesto | Vina | MM/GBSA total | ΔvdW + ΔEEL | ΔEGB | RMSD ligando |
|---|---|---|---|---|---|
| Quinine (neutra) | −7.69 | −55.59 | −60.63 | **+10.22** | **6.4 Å** |
| Pallidine | −8.07 | −42.12 | −56.05 | +18.80 | 5.5 Å |
| **Dabrafenib** | **−8.82** | −39.01 | **−75.82** | **+43.46** | 3.3 Å |
| Notholaenic acid | −7.72 | −36.67 | −55.91 | +24.68 | 2.6 Å |
| Macaridine | −7.38 | −32.54 | −47.08 | +19.14 | 2.8 Å |

Tres cosas se leen ahí:

1. **Por interacción directa (ΔvdW + ΔEEL) dabrafenib gana por 15 kcal/mol.** La física está bien capturada;
   lo que lo hunde es `EGB`, el término de desolvatación, que es la parte menos fiable del método. Conviene
   **reportar esa columna junto al total**.
2. **El primer puesto de la quinina es un artefacto.** Gana solo por pagar la desolvatación más baja, está
   modelada con la protonación equivocada, y su ligando se fue **6.4 Å del sitio acoplado** antes de empezar
   la producción: no estaba en el bolsillo de ATP cuando se calculó su energía.
3. **Pallidine también deriva** (5.5 Å) y sus puentes de hidrógeno caen del 94 % al 18 %.

**Faltan controles.** Ninguno de los cuatro productos naturales tiene actividad publicada frente a quinasas
(Cioffi *et al.* caracterizan el ácido notoloénico solo como antioxidante). Con un único control positivo se
puede demostrar que el método funciona técnicamente, pero no que **ordena** bien. Añadir vemurafenib como
segundo positivo y un inactivo como negativo es el mejor retorno por hora de GPU de la lista pendiente.

**10 ns no equivalen al paper base.** Singh *et al.* corrieron **900 ns**, con un RMSD medio de 3.05 Å en los
últimos 100 ns y tres conformaciones por agrupamiento. Aquí las pendientes de RMSD son ≤ 0.05 Å/ns, pero eso
es una meseta **local**: el sistema no ha hecho esas transiciones. Replicar 900 ns × 5 son 17 días de cola
en el mejor trozo de A100 y 36 en el peor. Más razonable: **100 ns con tres réplicas de semilla distinta**
(~2 días por réplica), que además mejora al paper base en un aspecto — ellos usaron una sola trayectoria por
compuesto, y sin réplicas no hay barras de error con las que defender que −55.6 y −42.1 son distintos.

**El orden correcto de trabajo** es: (1) protonaciones, (2) rehacer sus sistemas, (3) 10 ns de validación
—¿se queda la quinina protonada en el bolsillo?—, y solo entonces (4) decidir corridas largas. Lanzar
cientos de horas de GPU antes del paso 3 es gastarlas sobre moléculas con la carga equivocada.

**Protocolo frente al paper base.** Coincide en todo salvo la duración: ff19SB, GAFF/Antechamber, TIP3P,
0.15 M NaCl, calentamiento 0→310 K en 500 ps, paso de 2 fs, Langevin con colisión 1, PME con condiciones
periódicas, Vina con *exhaustiveness* 100 y MM/GBSA con `MMPBSA.py` despreciando la entropía.

---

## Notas

- Versiones usadas en el reporte: AmberTools 26, OpenMM 8.6.1 (8.2 en el HPC), AutoDock Vina 1.2.7,
  RDKit 2026.03, ParmEd 4.3.1, PDB2PQR 3.6–3.8 y MDAnalysis 2.10–2.11.
- **OpenMM 8.6 no sirve en CEDIA**: su driver es CUDA 12.2 y esos builds piden ≥ 12.9; fallan con
  *unsupported PTX*. Ahí va `environment-cuda12.yml` (OpenMM 8.2). En el portátil, con driver 595 / CUDA
  13.2, 8.6.1 funciona.
- Los nombres de carpeta **no llevan espacios** a propósito: las directivas `#SBATCH` se parten por espacios
  y no admiten comillas, así que `--output=.../10nsNotholaenic acid/...` rompe el envío.
- `resultados_50ns/` dejó de versionarse en el commit `d1c9729`. Los 154 archivos siguen en el historial y
  se recuperan con `git checkout d1c9729^ -- resultados_50ns && git reset -q resultados_50ns`.
- Los `correr_50ns.sh` de cada compuesto quedan como referencia histórica; lo que hacían lo cubre
  `dmlrun run --pipeline resultados`, que además elige el entorno según la PC y no fija el tiempo en el
  nombre de las carpetas.
