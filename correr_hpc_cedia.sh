#!/usr/bin/env bash
# =============================================================================
#  Lanza el pipeline de dinámica en la HPC de CEDIA (Slurm, perfil `hpc`).
#
#  Uso:
#      ./correr_hpc_cedia.sh [tiempo] [compuestos]
#      ENVIAR=1 ./correr_hpc_cedia.sh [tiempo] [compuestos]
#
#      tiempo      prueba | 10 | 20 | 50 | <ns>     (por defecto: 10)
#      compuestos  lista separada por comas         (por defecto: los cinco
#                                                    que tienen PDB en Data/)
#
#  Sin ENVIAR=1 solo comprueba el entorno y enseña el plan: NO encola nada.
#  Esa es la red de seguridad; revisa la salida antes de poner ENVIAR=1.
#
#  Antes de la primera ejecución:   git pull origin main
# =============================================================================
set -euo pipefail

TIEMPO="${1:-10}"
COMPUESTOS="${2:-dabrafenib,macaridine,pallidine,quinine,notholaenic_acid}"
SERVIDOR=hpc

RAIZ="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$RAIZ"

aviso() { printf '\n== %s\n' "$*"; }
malo()  { printf '   FALLA  %s\n' "$*"; ERRORES=$((ERRORES+1)); }
bien()  { printf '   ok     %s\n' "$*"; }
nota()  { printf '   aviso  %s\n' "$*"; }
ERRORES=0

# -----------------------------------------------------------------------------
#  Entorno. micromamba vive en ~/.local/bin y no siempre está en el PATH de un
#  shell no interactivo; sin esto el `activate` falla con "command not found".
# -----------------------------------------------------------------------------
aviso "Entorno"
export PATH="$HOME/.local/bin:$PATH"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"

if ! command -v micromamba >/dev/null; then
    printf '   FALLA  no encuentro el binario micromamba\n'
    printf '          búscalo con: find $HOME -maxdepth 3 -name micromamba -type f\n'
    exit 1
fi
eval "$(micromamba shell hook -s bash)"
micromamba activate md-braf
bien "entorno md-braf activo  ($(command -v python))"

# -----------------------------------------------------------------------------
#  OpenMM contra el driver. En CEDIA el driver es CUDA 12.2, y los builds de
#  OpenMM 8.6 exigen CUDA >= 12.9: fallan con "unsupported PTX". El entorno
#  correcto ahí es el de environment-cuda12.yml, que trae OpenMM 8.2.
# -----------------------------------------------------------------------------
aviso "OpenMM y GPU"
OPENMM="$(python -c 'import openmm; print(openmm.__version__)' 2>/dev/null || echo ausente)"
if [ "$OPENMM" = ausente ]; then malo "OpenMM no se puede importar"; else bien "OpenMM $OPENMM"; fi

if command -v nvidia-smi >/dev/null; then
    DRIVER="$(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1 || true)"
    if [ -n "$DRIVER" ]; then bien "driver NVIDIA $DRIVER"
    else printf '          (este nodo no expone GPU: normal en el de login)\n'; fi
fi
case "$OPENMM" in
    8.6*|8.7*|8.8*|8.9*|9.*)
        nota "OpenMM $OPENMM con driver CUDA 12.2 falla con PTX no soportado"
        printf '          usa el entorno de environment-cuda12.yml (OpenMM 8.2)\n';;
esac

# -----------------------------------------------------------------------------
#  Slurm
# -----------------------------------------------------------------------------
aviso "Slurm"
for exe in sbatch squeue scancel; do
    if command -v "$exe" >/dev/null; then bien "$exe"; else malo "falta $exe"; fi
done

# La QOS 'normal' de CEDIA permite UN trabajo corriendo a la vez. Una sesión de
# escritorio de Open OnDemand consume ese cupo y deja todo lo demás en espera.
EN_MARCHA="$(squeue -u "$USER" -h -t RUNNING -o '%i %j' 2>/dev/null || true)"
if [ -n "$EN_MARCHA" ]; then
    printf '   trabajos tuyos ya en marcha (ocupan tu cupo de la QOS normal):\n'
    printf '%s\n' "$EN_MARCHA" | sed 's/^/        /'
    if printf '%s' "$EN_MARCHA" | grep -qi 'dashboard\|desktop\|bc_'; then
        nota "una sesión de Open OnDemand ocupa tu único cupo"
        printf '          entra por Clusters -> Shell Access (no consume cola) y ciérrala\n'
        printf '          con scancel, o lo que encoles ahora esperará a que termine\n'
    fi
else
    bien "sin trabajos en marcha: hay cupo libre"
fi

# -----------------------------------------------------------------------------
#  Herramientas de la etapa `sistema` (docking + topología AMBER). Solo hacen
#  falta para compuestos sin sistema/sistema.prmtop; si todos vienen preparados
#  del repositorio, esa etapa se omite y da igual que falten.
# -----------------------------------------------------------------------------
aviso "Preparación del sistema"
FALTA_SISTEMA=0
IFS=',' read -ra LISTA <<< "$COMPUESTOS"
for c in "${LISTA[@]}"; do
    DIR="$(python - "$c" <<'PY' 2>/dev/null || true
import sys
from dmlrun.config import cargar_proyecto
print(cargar_proyecto().compuesto(sys.argv[1]).sistema_dir)
PY
)"
    if [ -n "$DIR" ] && [ -f "$DIR/sistema.prmtop" ]; then
        bien "$c: topología ya preparada (se omite el docking)"
    else
        printf '   %s: hay que preparar el sistema\n' "$c"
        FALTA_SISTEMA=1
    fi
done

if [ "$FALTA_SISTEMA" = 1 ]; then
    for exe in obabel vina antechamber parmchk2 tleap; do
        if command -v "$exe" >/dev/null; then bien "$exe"
        else malo "falta $exe (lo necesita la etapa sistema)"; fi
    done
fi

# -----------------------------------------------------------------------------
#  Plan
# -----------------------------------------------------------------------------
if [ "$ERRORES" -gt 0 ]; then
    printf '\n%s comprobación(es) fallida(s). No se envía nada.\n' "$ERRORES"
    exit 1
fi

aviso "Plan ($TIEMPO ns - $COMPUESTOS)"
dmlrun run --pipeline resultados --compuesto "$COMPUESTOS" \
           --tiempo "$TIEMPO" --servidor "$SERVIDOR" --dry-run

if [ "${ENVIAR:-0}" != 1 ]; then
    printf '\nNada enviado. Revisa el plan de arriba y, si está bien:\n\n'
    printf '    ENVIAR=1 %s %s %s\n\n' "$0" "$TIEMPO" "$COMPUESTOS"
    printf 'Con MaxJobs=1 los trabajos se ejecutan de uno en uno; Slurm los\n'
    printf 'encola y los va lanzando solo.\n'
    exit 0
fi

aviso "Enviando"
dmlrun run --pipeline resultados --compuesto "$COMPUESTOS" \
           --tiempo "$TIEMPO" --servidor "$SERVIDOR"

aviso "Cola"
squeue -u "$USER" -o "%.8i %.24j %.9T %.10M %.28R"

printf '\nSeguimiento:\n'
printf '    dmlrun status\n'
printf '    squeue -u $USER\n'
printf '    tail -f resultados/<sesion>/<tiempo>ns<Compuesto>/estado/ejecucion.log\n\n'
printf 'Los resultados van a resultados/<sesion>/<tiempo>ns<Compuesto>/ (md,\n'
printf 'resultados, analisis, estado): una carpeta por tanda de lanzamiento. Es\n'
printf 'reanudable: si Slurm corta por tiempo, relanza y sigue donde iba.\n'
