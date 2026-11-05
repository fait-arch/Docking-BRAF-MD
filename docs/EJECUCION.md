# Ejecución de los pipelines en HPC y en otros servidores Linux

Toda la orquestación se hace con una sola CLI, `dmlrun`. Lo único que cambia
entre entornos es `--servidor`: el gestor de colas, los recursos, los módulos,
el entorno de conda y la plataforma de OpenMM salen de
[`config/servidores.yaml`](../config/servidores.yaml). No hay que editar código
ni los scripts de simulación para cambiar de máquina ni de tiempo de simulación.

```
dmlrun servidores     perfiles de entorno disponibles
dmlrun compuestos     compuestos registrados y sus ejecuciones
dmlrun config         muestra y comprueba un perfil en la máquina actual
dmlrun run            envía ejecuciones
dmlrun status         estado y progreso
dmlrun logs           log de una ejecución
dmlrun cancel         cancela una ejecución
dmlrun comparar       compara tiempos de simulación de un compuesto
```

---

## 1. Instalación en un servidor nuevo

Los mismos cuatro pasos sirven para la HPC y para cualquier otro Linux.

```bash
git clone <url-del-repo> dockingML
cd dockingML

# entorno (AmberTools, OpenMM, MDAnalysis, RDKit, Vina, PyYAML…)
conda env create -f simulaciones/Pallidine_24543347/environment.yml
conda activate md-braf

# la CLI y los comandos dml-* del pipeline de ML
pip install -e .

# comprueba que el perfil del servidor encaja con esta máquina
dmlrun config --servidor hpc          # en la HPC
dmlrun config --servidor linux-gpu    # en otro servidor con GPU
```

`pip install -e .` deja disponibles `dmlrun` y los comandos `dml-dock`,
`dml-features`, `dml-rescore` y `dml-train` que usa el pipeline de ML. Si no
quieres instalar nada, `python -m dmlrun.cli …` funciona igual para todo salvo
las etapas de ML.

Si ejecutas desde una carpeta que no es la raíz del repo:

```bash
export DMLRUN_RAIZ=/ruta/a/dockingML
```

---

## 2. Perfiles de servidor

```bash
dmlrun servidores
```

| perfil | cola | para qué |
|---|---|---|
| `hpc` | Slurm, `interactive`, 1 GPU | dinámica molecular en la HPC |
| `hpc-cpu` | Slurm, `interactive`, sin GPU | análisis, MM/GBSA y ML en la HPC |
| `linux-gpu` | directa (`setsid`) | servidor externo con GPU NVIDIA |
| `linux-cpu` | directa (`setsid`) | servidor externo sólo CPU |
| `local` | directa | máquina de trabajo, pruebas cortas |
| `actual` | directa | aquí mismo, con el entorno ya activo |

`--servidor auto` elige solo: respeta `$DMLRUN_SERVIDOR`, luego un perfil cuyo
nombre coincida con el hostname, luego el primer perfil de Slurm si existe
`sbatch`, y en último lugar `local`.

**Añadir un servidor** es copiar un bloque en `config/servidores.yaml`:

```yaml
  mi-servidor:
    descripcion: "Servidor del laboratorio"
    scheduler: local            # local | slurm
    cpus: 16
    memoria_gb: 64
    gpus: 1
    plataforma_openmm: CUDA
    conda_env: md-braf
    pre_comandos:
      - 'source "$HOME/miniconda3/etc/profile.d/conda.sh"'
    variables:
      CUDA_VISIBLE_DEVICES: "0"
```

Desde ese momento `--servidor mi-servidor` ya funciona.

---

## 3. Los dos pipelines

### `--pipeline resultados` (MD + análisis + MM/GBSA)

Etapas: `sistema`, `min`, `cal`, `eq`, `prod`, `rmsd`, `analisis`, `mmgbsa`.

### `--pipeline ml` (docking + descriptores + rescoring)

Etapas: `centro`, `docking`, `complejos`, `features`, `rescore`.

Con `--fuente docking` acopla de nuevo el ligando y puntúa esas poses: es la
referencia común, independiente del tiempo de simulación. Con `--fuente md`
toma fotogramas repartidos por la producción del tiempo que le indiques
(`--frames-ml 20`), de modo que se puede comparar el rescoring a 10, 20 y 50 ns
**sobre los mismos compuestos**. En ese caso la etapa `docking` no se usa.

`--pipeline ambos` lanza los dos, uno detrás de otro.

---

## 4. Comandos de ejecución

Todos aceptan `--dry-run` (muestra el plan sin enviar nada) y `--solo-script`
(escribe el `sbatch` en `runs/<tiempo>/estado/` sin enviarlo, para revisarlo o
enviarlo a mano).

### En la HPC (Slurm)

```bash
# los dos compuestos, a 10, 20 y 50 ns, en la cola con GPU
dmlrun run --pipeline resultados --compuesto todos --tiempo 10,20,50 --servidor hpc

# ver primero qué se enviaría y con qué parámetros
dmlrun run --pipeline resultados --compuesto todos --tiempo 10,20,50 \
       --servidor hpc --dry-run

# ver el script de Slurm completo
dmlrun run --pipeline resultados --compuesto pallidine --tiempo 50 \
       --servidor hpc --dry-run --mostrar-script

# sólo la producción, como máximo 2 segmentos por envío (trabajos más cortos)
dmlrun run --pipeline resultados --compuesto pallidine --tiempo 50 \
       --etapas prod --max-segmentos 2 --servidor hpc

# análisis y MM/GBSA sin gastar GPU
dmlrun run --pipeline resultados --compuesto todos --tiempo 50 \
       --etapas rmsd,analisis,mmgbsa --servidor hpc-cpu

# pipeline de ML sobre la trayectoria de 50 ns
dmlrun run --pipeline ml --compuesto todos --tiempo 50 \
       --fuente md --frames-ml 20 --servidor hpc-cpu

# ML de referencia (redocking, sin depender de la MD)
dmlrun run --pipeline ml --compuesto todos --tiempo 50 \
       --fuente docking --motor vina --servidor hpc-cpu

# con un modelo entrenado
dmlrun run --pipeline ml --compuesto todos --tiempo 20 \
       --fuente md --modelo modelos/rescorer_xgb.pkl --servidor hpc-cpu
```

Recursos del perfil sobrescritos sólo para este envío:

```bash
dmlrun run --pipeline resultados --compuesto pallidine --tiempo 50 \
       --servidor hpc --cpus 16 --memoria 48 --horas 24
```

### En otro servidor Linux (sin cola)

Exactamente los mismos comandos, cambiando el perfil. El trabajo se lanza con
`setsid`, así que sobrevive al cierre de la terminal.

```bash
dmlrun run --pipeline resultados --compuesto todos --tiempo 20 --servidor linux-gpu
dmlrun run --pipeline ml --compuesto pallidine --tiempo 20 --fuente md --servidor linux-cpu
```

### Prueba rápida antes de ocupar la cola

```bash
dmlrun run --pipeline resultados --compuesto pallidine --tiempo prueba --servidor local
```

`prueba` son 100 ps con el mismo protocolo (20 ps de calentamiento, 20 de
equilibrado, 4 segmentos y un fotograma cada ps).

---

## 5. Seguimiento

```bash
dmlrun status                                  # últimas 15 ejecuciones
dmlrun status --todas
dmlrun status --compuesto pallidine
dmlrun status --id pallidine-resultados-50ns   # un prefijo del id basta
dmlrun logs   --id pallidine-resultados-50ns -n 100
dmlrun logs   --id pallidine-resultados-50ns -f     # en vivo
dmlrun cancel --id pallidine-resultados-50ns
```

La columna **progreso** no viene del gestor de colas: se lee del disco
(`producción 7/10` cuenta los segmentos con su `.xml` guardado), así que es
correcta aunque Slurm haya matado el trabajo por tiempo.

El estado funciona igual con Slurm y sin cola porque cada trabajo escribe su
propio código de salida en `runs/<tiempo>/estado/codigo_salida`.

---

## 6. Dónde quedan los resultados

```
simulaciones/<Compuesto>/
├── sistema/                     topología AMBER: se prepara una vez y vale
│                                para todos los tiempos de simulación
└── runs/
    ├── 10ns/
    │   ├── md/                  01_min, 02_cal, 03_eq, 04_prod_00N
    │   ├── analisis/            termodinámica, RMSD, contactos, mmgbsa/
    │   ├── resultados/          figuras
    │   ├── ml/                  complejos, features.csv, rescored.csv
    │   └── estado/              log, script, manifest.json, codigo_salida
    ├── 20ns/   …
    └── 50ns/   …
```

**Compatibilidad con lo que ya existe**: las carpetas antiguas
`md_prueba`/`analisis`/`resultados_prueba` y `md_50ns`/`analisis_50ns`/
`resultados_50ns` se siguen usando tal cual; no hay que mover ni duplicar nada.
`dmlrun` las detecta y las reutiliza (`(carpetas heredadas)` en la salida), y
también reconoce las etapas ya completadas dentro de ellas. Fuerza una u otra
estructura con `--layout runs` o `--layout legacy`.

---

## 7. Sacar los resultados: `dmlrun reporte`

Las figuras y las tablas de una corrida quedan repartidas entre varias carpetas.
Este comando las reune en una sola, por compuesto y tiempo, y escribe un informe
que puede leer alguien que no lanzo la simulacion:

```bash
dmlrun reporte                                   # todas las corridas que existan
dmlrun reporte --compuesto pallidine --tiempo 50
dmlrun reporte --destino ~/informes_braf         # otra carpeta de salida
```

Resultado:

```
reportes/
├── indice.html                  tabla comparativa con enlaces a cada informe
└── pallidine_50ns/
    ├── informe.html             autocontenido: imagenes incrustadas
    ├── figuras/                 PNG y PDF vectorial de cada figura
    └── datos/                   CSV, resumen.json, FINAL_MMGBSA.dat, manifest
```

`informe.html` lleva las imagenes dentro del propio archivo, asi que se abre sin
conexion, se copia a un pendrive o se envia por correo y se ve igual en
cualquier ordenador. Tambien se imprime a PDF desde el navegador.

Cada figura va acompanada de **que muestra** y **que hay que mirar en ella**
(que el RMSD de Ca llegue a una meseta, que residuos superan el 90 % de
ocupacion y son candidatos a mutagenesis, por que el MM/GBSA no es una afinidad
absoluta al faltarle el termino entropico...). Los numeros van en tablas con sus
unidades y una columna de lectura.

El informe avisa en su cabecera cuando la corrida es demasiado corta para
interpretarse, y el indice marca esas filas: una prueba de 100 ps sirve para
comprobar el montaje, no para sacar conclusiones biologicas.

Los `figuras/*.pdf` son vectoriales y valen para publicacion; los `datos/*.csv`
permiten rehacer o reestilizar cualquier grafica.

---

## 8. Reproducibilidad entre entornos

Cada envío escribe `runs/<tiempo>/estado/manifest.json` con:

- **huella**: hash de los parámetros científicos (compuesto, SMILES, carga, pH,
  residuos, tiempo, segmentos, semilla…). **Misma huella = resultados
  comparables**, aunque se hayan calculado en servidores distintos. Cambiar el
  tiempo de simulación cambia la huella a propósito;
- commit de git, rama y si el árbol estaba sucio;
- versiones de OpenMM, ParmEd, MDAnalysis, mdtraj, NumPy, RDKit, scikit-learn…;
- rutas de `tleap`, `antechamber`, `MMPBSA.py`, `vina`, `obabel`;
- servidor, scheduler, plataforma de OpenMM, CPU/GPU y hostname.

La semilla del integrador es fija (`--semilla 2022`) y entra en la huella, así
que dos ejecuciones del mismo compuesto y tiempo en la misma plataforma son
reproducibles. Entre plataformas distintas (CUDA frente a CPU) OpenMM no
garantiza trayectorias bit a bit idénticas; el manifiesto deja constancia de la
plataforma usada para que la diferencia sea explicable y no una sorpresa.

Comparar lo que hay:

```bash
dmlrun comparar --compuesto pallidine
dmlrun comparar --compuesto todos --tiempo 10,20,50
```

```
  compuesto   tiempo  producción  mmgbsa  ml   huella        servidor
  ----------  ------  ----------  ------  ---  ------------  --------
  pallidine   50ns    8/10        sí      sí   0653f837ae79  hpc
  pallidine   20ns    4/4         sí      no   a91c40bd1e72  linux-gpu
```

---

## 9. Añadir un compuesto

1. Pon sus PDB de receptor y ligando en `Simulacion/Data`.
2. Añade un bloque en `config/compuestos.yaml`:

```yaml
  mi_compuesto:
    nombre: MiCompuesto
    pubchem: 12345678
    directorio: simulaciones/MiCompuesto_12345678
    titulo: "MiCompuesto (12345678) – BRAF V600E"
    receptor: MiCompuesto12345678_receptor.pdb
    ligando: MiCompuesto12345678_ligand.pdb
    smiles: "..."
    carga_ligando: 0
```

3. Prepara el sistema y lanza:

```bash
dmlrun run --pipeline resultados --compuesto mi_compuesto --tiempo prueba \
       --etapas sistema --servidor local
dmlrun run --pipeline resultados --compuesto mi_compuesto --tiempo 10,20,50 \
       --servidor hpc
```

No hay que copiar scripts: los de `simulaciones/scripts/` son los mismos para
todos los compuestos.

---

## 10. Reanudación y problemas frecuentes

**Reanudar.** Los envíos son idempotentes: `dmlrun run` omite las etapas que ya
tienen su salida (`min`, `cal`, `eq`, cada segmento de producción, `mmgbsa`).
Si un trabajo se corta por límite de tiempo, vuelve a lanzar el mismo comando y
continúa donde estaba. Con `--forzar` repite todo.

**`sbatch: Requested node configuration is not available`.** Ningún nodo de la
partición existe con esa combinación. `dmlrun` lo explica al fallar; revisa en
este orden:

```bash
scontrol show node compute-0-0 | grep -i gres    # nombre real del GRES de GPU
sinfo -p interactive -o "%P %c %m %G"            # CPU, memoria y GPU por nodo
scontrol show partition interactive              # límite de tiempo
```

Ajusta `gres`, `cpus`, `memoria_gb` u `horas` en el perfil. En nodos con GPU
particionada (MIG) el GRES no es `gpu:1` sino algo como `gpu:2g.10gb`.

**`conda: command not found` en el nodo de cómputo.** El nodo no hereda tu
`.bashrc`. Pon la ruta real en `pre_comandos` del perfil.

**MM/GBSA.** El número de fotogramas y el intervalo se calculan contando los
fotogramas que **realmente** hay en los `.dcd`, no el tiempo nominal: funciona
igual si la producción quedó a medias.
