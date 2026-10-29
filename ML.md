# Versión con machine learning (ML/MM con ANI-2x)

Esta versión repite **el mismo protocolo de dinámica molecular** que la versión clásica (sección 2.1.3 de
*RSC Adv.* 2022, 12, 30181) con un solo cambio: la energía interna del ligando se calcula con el potencial
de red neuronal **ANI-2x** en lugar de con GAFF. Sirve para comparar, en el mismo sistema y con la misma
duración, qué cambia en el comportamiento del ligando cuando su descripción viene de un modelo de machine
learning entrenado con cálculos cuánticos (DFT).

## Qué se describe con ML y qué no

| Parte del sistema | Versión clásica (`correr_50ns.sh`) | Versión ML/MM (`correr_ml.sh`) |
|---|---|---|
| Enlaces, ángulos, torsiones e interacciones internas del ligando | GAFF + cargas AM1-BCC | **ANI-2x** (conjunto de 8 redes neuronales) |
| Interacciones ligando–proteína y ligando–agua | Lennard-Jones + cargas (GAFF/AM1-BCC) | Igual (acoplamiento mecánico) |
| Proteína, agua, iones | ff19SB, TIP3P | Igual |
| Caja, sal, 310 K, 1 atm, Langevin γ = 1 ps⁻¹, 2 fs, corte 11 Å, PME | según el artículo | Igual |
| Minimización (10 000 it., k = 10), calentamiento 500 ps, equilibrado 500 ps | según el artículo | Igual, con el ligando en ANI-2x |
| Enlaces con H del ligando | rígidos | rígidos (para conservar el paso de 2 fs) |

El acoplamiento es **mecánico**: la red neuronal da la geometría y la flexibilidad del ligando (torsiones,
anillos, conjugación), pero la interacción con la proteína sigue usando las cargas y el Lennard-Jones del
campo de fuerza. Es el esquema estándar de `openmm-ml` (Eastman *et al.*, *J. Phys. Chem. B* 2024).

**Limitación importante:** ANI-2x solo está entrenado con moléculas **neutras** con H, C, N, O, F, S y Cl.
Sirve para Macaridine y para Pallidine neutra; **no** para Pallidine protonada (NH⁺, carga +1). El script se
detiene con un mensaje si el ligando tiene carga.

## Instalación (una vez por PC con GPU NVIDIA)

```bash
cd Docking-BRAF-MD
git pull
micromamba create -f environment-ml-cuda12.yml   # entorno md-braf-ml recomendado (con NNPOps); driver CUDA >= 12.2
# alternativa con las versiones más nuevas, pero sin NNPOps (unas 5 veces más lenta en GPU):
#   micromamba create -f environment-ml.yml
```

**NNPOps** reemplaza las operaciones de TorchANI por núcleos CUDA optimizados para la misma red ANI-2x (mismas
energías y fuerzas). El artículo de OpenMM 8 midió que es unas 5.7 veces más rápido. Solo existe para
PyTorch 2.5 + TorchANI 2.2, por eso `environment-ml-cuda12.yml` fija esas versiones. Los scripts lo usan
automáticamente si está instalado (`IMPL=auto`). Si el entorno se creó antes sin NNPOps:

```bash
micromamba install -y -n md-braf-ml -c conda-forge "nnpops=0.6" "pytorch=2.5" "torchani=2.2.4" \
                   "openmm=8.2" "openmm-torch=1.5" "cuda-version=12.2"
```

`MAMBA_ROOT_PREFIX` debe ser una ruta **absoluta** (por ejemplo `$HOME/micromamba`). Si vale algo relativo,
como `y`, micromamba busca los entornos en una carpeta distinta según el directorio actual y `correr_ml.sh`
se detiene con un aviso.

El entorno `md-braf` (el de siempre) se sigue usando para el análisis y el MM/GBSA; `correr_ml.sh` llama a
cada entorno cuando lo necesita, así que no hace falta activar ninguno.

Comprobación rápida (construye el sistema ML/MM, evalúa energía y fuerzas y dice si ve la GPU):

```bash
cd simulaciones/Macaridine_94717379/scripts
micromamba run -n md-braf-ml python potencial_ml.py --sistema ../sistema
```

Debe terminar con `Verificación correcta` y `CUDA disponible: True`, e indica la implementación usada y el
costo de la red en milisegundos por paso. La primera vez descarga los pesos de ANI-2x de GitHub.

**Elegir la variante antes de lanzar** (unos minutos; compara las cuatro combinaciones en la GPU):

```bash
for imp in nnpops torchani; do for red in -1 0; do echo "== $imp, red $red"
  micromamba run -n md-braf-ml python potencial_ml.py --sistema ../sistema \
      --ml-implementacion $imp --ml-modelo-indice $red 2>&1 | grep -E "Costo|ERROR"; done; done
```

Con c milisegundos por paso de la red y ~1.4 ms del resto del sistema (A100 MIG de 10 GB), la velocidad es
≈ 172.8 / (1.4 + c) ns/día. Para que 20 ns quepan en una sesión de 16 h hace falta c ≲ 4.5 ms; para
dos sesiones, c ≲ 10 ms. Medido en el HPC de CEDIA con TorchANI y el conjunto de 8 redes: c ≈ 41 ms
(4 ns/día), demasiado lento.

## Correr: prueba de 20 ns y después 50 ns

```bash
cd simulaciones/Macaridine_94717379/scripts          # o Pallidine_24543347 (forma neutra)

# 1) prueba de 20 ns (minimización + 0.5 ns + 0.5 ns + 20 ns de producción, análisis, MM/GBSA y comparación)
nohup bash correr_ml.sh > corrida_ml_20ns.log 2>&1 &

# 2) cuando termine, ampliar la MISMA corrida a 50 ns: reutiliza los 20 ns y corre solo los 30 ns que faltan
T_PROD_PS=50000 nohup bash correr_ml.sh > corrida_ml_50ns.log 2>&1 &
```

Ver el avance: `tail -f corrida_ml_20ns.log` (o `tail -f ../md_ml_ani2x/run.log`).

Antes de la producción, el script mide la velocidad y la guarda en `../md_ml_ani2x/velocidad_auto.txt`, con la
estimación de horas para 20 y 50 ns. Como referencia, el artículo de OpenMM 8 midió que describir un ligando
con ANI-2x hace la simulación entre 2.7 y ~30 veces más lenta que la clásica, según la implementación. Si va
demasiado lenta, hay dos opciones (cada una crea su propia carpeta, para no mezclar corridas):

```bash
RED=0 bash correr_ml.sh                       # una sola red en vez del conjunto de 8: ~2-3 veces más rápido
EXTRA_ML="--ml-cuda-graphs" bash correr_ml.sh # CUDA graphs en la GPU (si la versión de openmm-torch lo admite)
```

## Qué produce

| Carpeta (en `simulaciones/<Compuesto>/`) | Contenido |
|---|---|
| `md_ml_ani2x/` | trayectorias (`.dcd`), registros, estados (`.xml`), `metodo.txt`, `velocidad_auto.txt` |
| `analisis_ml_ani2x_20ns/` y `_50ns/` | mismos análisis que la versión clásica (RMSD, RMSF, contactos, puentes de H, MM/GBSA) |
| `resultados_ml_ani2x_20ns/` y `_50ns/` | mismas figuras que la versión clásica |
| `comparacion_ml_ani2x_20ns/` y `_50ns/` | **comparación clásico frente a ML/MM** en la misma ventana de tiempo |

La comparación (`06_comparar_clasico_ml.py`) necesita la corrida clásica en `../md_50ns` (y su MM/GBSA en
`../analisis_50ns/mmgbsa`). Si no está, primero `bash correr_50ns.sh`. Calcula para las dos corridas:
RMSD de Cα y del ligando, radio de giro, puentes de H (ocupación y número medio), contactos por residuo,
**torsiones del ligando** (distribuciones y su solapamiento, que es lo que más cambia de GAFF a ANI-2x) y
MM/GBSA. Escribe `tabla_comparacion.csv`, `torsiones_ligando.csv`, `resumen_comparacion.json` y cuatro
figuras (`fig_cmp_*.png/pdf`).

## Cómo interpretar la comparación

- **Una réplica por método no basta para afirmar diferencias.** El momento en que Macaridine pierde el puente
  con Thr529 (a los ~28 ns en la corrida clásica) puede cambiar de una corrida a otra aunque el método sea el
  mismo. Para una conclusión firme, correr al menos 3 réplicas de cada versión (cambiando `--semilla`).
- **El MM/GBSA se calcula igual en las dos** (con el campo de fuerza): sus diferencias reflejan las poses y
  conformaciones que visita cada simulación, no la energía de ANI-2x.
- **La energía potencial de los registros no es comparable** entre versiones: en la ML/MM se cambia la
  energía interna del ligando (se resta la de ANI-2x en la geometría inicial para mantener la escala).
- Lo que la versión ML/MM puede mejorar es la **geometría interna del ligando**. En Pallidine, dos diedros del
  metoxilo de la dienona tenían parámetros genéricos de GAFF (penalización 237); ahí es donde más se espera
  ver diferencias.

## Archivos de la versión ML

- `environment-ml.yml`, `environment-ml-cuda12.yml`: entorno `md-braf-ml`.
- `scripts/potencial_ml.py`: convierte el sistema clásico en ML/MM (TorchForce con ANI-2x) y lo verifica.
- `scripts/02_dinamica.py --ml ani2x`: misma dinámica, con el ligando en ANI-2x.
- `scripts/benchmark.py --ml ani2x`: velocidad de la versión ML/MM.
- `scripts/correr_ml.sh`: todo el flujo (dinámica, análisis, MM/GBSA y comparación).
- `scripts/06_comparar_clasico_ml.py`: comparación clásico frente a ML/MM.
