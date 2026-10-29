#!/usr/bin/env bash
# =============================================================================
#  Versión con machine learning (ML/MM) del protocolo del artículo (RSC Adv. 2022, 12, 30181).
#  Todo es igual que en correr_50ns.sh (ff19SB, TIP3P, caja, NaCl 0.15 M, 310 K, 1 atm, Langevin,
#  2 fs, corte 11 Å, PME, minimización, calentamiento y equilibrado de 500 ps), salvo que la energía
#  interna del ligando se calcula con el potencial de red neuronal ANI-2x en lugar de GAFF.
#
#  Uso (desde esta carpeta; no hace falta activar ningún entorno):
#     bash correr_ml.sh                       # prueba de 20 ns
#     T_PROD_PS=50000 bash correr_ml.sh       # amplía la MISMA corrida a 50 ns (reutiliza los 20 ns)
#  En segundo plano:
#     nohup bash correr_ml.sh > corrida_ml_20ns.log 2>&1 &
#
#  Necesita dos entornos de micromamba:
#     md-braf-ml  dinámica ML/MM (environment-ml.yml de la raíz del repositorio)
#     md-braf     análisis y MM/GBSA (environment.yml)
#  Es reanudable. Al final compara con la corrida clásica de ../md_50ns en la misma ventana de tiempo.
#  Solo para ligandos neutros: ANI-2x no admite moléculas con carga.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

T_CAL_PS=${T_CAL_PS:-500}           # calentamiento (ps), igual que el artículo
T_EQ_PS=${T_EQ_PS:-500}             # equilibrado (ps), igual que el artículo
T_PROD_PS=${T_PROD_PS:-20000}       # producción: 20 ns (prueba); 50000 para 50 ns
SEG_PS=${SEG_PS:-5000}              # segmentos de 5 ns, como la versión clásica
GUARDAR_PS=${GUARDAR_PS:-10}        # un fotograma cada 10 ps
MMGBSA_CADA_PS=${MMGBSA_CADA_PS:-250}   # MM/GBSA con un fotograma cada 250 ps, como la clásica
MODELO=${MODELO:-ani2x}             # ani2x (H C N O F S Cl) o ani1ccx (H C N O)
RED=${RED:--1}                      # -1 = conjunto de 8 redes; 0..7 = una sola red (~2-3 veces más rápido)
IMPL=${IMPL:-auto}                  # auto = NNPOps si está instalado (mucho más rápido en GPU); torchani; nnpops
MIN_ITER=${MIN_ITER:-10000}
BENCH_PASOS=${BENCH_PASOS:-1000}    # 0 = no medir la velocidad
EXTRA_ML=${EXTRA_ML:-}              # opciones extra de 02_dinamica.py, p. ej. "--ml-cuda-graphs"
ENV_ML=${ENV_ML:-md-braf-ml}
ENV_MD=${ENV_MD:-md-braf}
CLASICO=${CLASICO:-../md_50ns}
GB_CLASICO=${GB_CLASICO:-../analisis_50ns/mmgbsa}

MM=${MAMBA_EXE:-$(command -v micromamba || true)}
[ -n "$MM" ] || { echo "ERROR: no encuentro micromamba (abrir una terminal nueva o definir MAMBA_EXE)"; exit 1; }
case "${MAMBA_ROOT_PREFIX:-}" in
    /*) ;;
    "") echo "AVISO: MAMBA_ROOT_PREFIX no está definido; micromamba usará su carpeta por defecto." ;;
    *)  echo "ERROR: MAMBA_ROOT_PREFIX='$MAMBA_ROOT_PREFIX' es una ruta relativa: micromamba buscaría los entornos"
        echo "       en una carpeta distinta según el directorio actual. Definirlo con una ruta absoluta, p. ej."
        echo "       export MAMBA_ROOT_PREFIX=\"\$HOME/micromamba\"  (y lo mismo en ~/.bashrc)"; exit 1 ;;
esac
en_ml() { "$MM" run -n "$ENV_ML" "$@"; }
en_md() { "$MM" run -n "$ENV_MD" "$@"; }
for e in "$ENV_ML" "$ENV_MD"; do
    "$MM" run -n "$e" true 2>/dev/null || { echo "ERROR: no existe el entorno '$e' en ${MAMBA_ROOT_PREFIX:-la carpeta por defecto} (ver ML.md)"; exit 1; }
done
calc() { awk "BEGIN{print $1}"; }

COMPUESTO=$(basename "$(cd .. && pwd)")
SEGMENTOS=$(awk -v t="$T_PROD_PS" -v s="$SEG_PS" 'BEGIN{n=int(t/s+0.5); if(n<1)n=1; print n}')
T_NS=$(awk -v t="$T_PROD_PS" 'BEGIN{printf "%g", t/1000}')
SUFIJO=$MODELO; [ "$RED" -ge 0 ] && SUFIJO=${MODELO}_red${RED}
MD=../md_ml_${SUFIJO}
ANA=../analisis_ml_${SUFIJO}_${T_NS}ns
FIG=../resultados_ml_${SUFIJO}_${T_NS}ns
CMP=../comparacion_ml_${SUFIJO}_${T_NS}ns
mkdir -p "$MD" "$ANA" "$FIG"

# la carpeta de trayectoria debe usar siempre la misma longitud de segmento (para poder ampliarla)
if [ -f "$MD/segmento_ps.txt" ] && [ "$(cat "$MD/segmento_ps.txt")" != "$SEG_PS" ]; then
    echo "ERROR: $MD usa segmentos de $(cat "$MD/segmento_ps.txt") ps y se pidió SEG_PS=$SEG_PS"; exit 1
fi
echo "$SEG_PS" > "$MD/segmento_ps.txt"
ML_OPC="--ml $MODELO --ml-modelo-indice $RED --ml-implementacion $IMPL $EXTRA_ML"

echo "== $COMPUESTO · ML/MM con $MODELO · producción $T_NS ns ($SEGMENTOS segmentos de $SEG_PS ps) -> $MD"
echo "== 0) Comprobación del potencial ML (carga del ligando, elementos, GPU)"
en_ml python potencial_ml.py --sistema ../sistema --ml "$MODELO" --ml-modelo-indice "$RED" \
      --ml-implementacion "$IMPL" --out "$MD/verificacion"

ETAPAS="min,cal,eq,prod"
[ -f "$MD/01_min.xml" ] && ETAPAS=${ETAPAS/min,/}
[ -f "$MD/02_cal.xml" ] && ETAPAS=${ETAPAS/cal,/}
[ -f "$MD/03_eq.xml" ]  && ETAPAS=${ETAPAS/eq,/}
DIN="python 02_dinamica.py $ML_OPC --sistema ../sistema --out $MD --min-iter $MIN_ITER \
     --t-calentamiento $T_CAL_PS --t-equilibrado $T_EQ_PS --t-produccion $T_PROD_PS \
     --segmentos $SEGMENTOS --guardar-cada $GUARDAR_PS"

if [[ "$ETAPAS" == min,* ]]; then
    echo "== 1) Minimización ML/MM"
    en_ml $DIN --etapas min | tee -a "$MD/run.log"
    ETAPAS=${ETAPAS/min,/}
fi
if [ "$BENCH_PASOS" -gt 0 ] && [ ! -f "$MD/velocidad_${IMPL}.txt" ]; then
    echo "== 2) Velocidad esperada (ML/MM)"
    en_ml python benchmark.py --sistema ../sistema --estado "$MD/01_min.xml" --pasos "$BENCH_PASOS" $ML_OPC \
        | tee "$MD/velocidad_${IMPL}.txt" || true
fi

echo "== 3) Dinámica ML/MM ($ETAPAS)"
en_ml $DIN --etapas "$ETAPAS" | tee -a "$MD/run.log"

N_HECHOS=$( (ls "$MD"/04_prod_*.xml 2>/dev/null || true) | wc -l)
if [ "$N_HECHOS" -gt "$SEGMENTOS" ]; then
    echo "AVISO: $MD tiene $N_HECHOS segmentos (más de los $SEGMENTOS pedidos); el análisis los incluye todos."
fi

echo "== 4) RMSD y análisis completo"
en_md python 03_rmsd.py --sistema ../sistema --md "$MD" --out "$FIG" \
      --titulo "$COMPUESTO – BRAF V600E · ML/MM ($MODELO)"
en_md python 04_analisis_reporte.py --base .. --md "$(basename "$MD")" --fig "$FIG" --out "$ANA"

echo "== 5) MM/GBSA (mismo modelo GB que la versión clásica, un fotograma cada $MMGBSA_CADA_PS ps)"
mkdir -p "$ANA/mmgbsa"
N_FRAMES=$(awk -v t="$T_PROD_PS" -v g="$GUARDAR_PS" 'BEGIN{print int(t/g+0.5)}')
INTERVALO=$(awk -v c="$MMGBSA_CADA_PS" -v g="$GUARDAR_PS" 'BEGIN{n=int(c/g+0.5); if(n<1)n=1; print n}')
TOPS=""
for d in ../analisis_50ns/mmgbsa ../analisis/mmgbsa; do
    [ -f "$d/com.prmtop" ] && [ -f "$d/rec.prmtop" ] && [ -f "$d/lig.prmtop" ] && { TOPS=$d; break; }
done
if [ -n "$TOPS" ]; then
    cp "$TOPS/com.prmtop" "$TOPS/rec.prmtop" "$TOPS/lig.prmtop" "$ANA/mmgbsa/"
else
    ( cd "$ANA/mmgbsa" && rm -f com.prmtop rec.prmtop lig.prmtop \
      && printf 'changeRadii mbondi2\noutparm com.prmtop\nquit\n' \
         | en_md parmed -n -p ../../sistema/complejo_seco.prmtop > parmed.log 2>&1 \
      && en_md ante-MMPBSA.py -p com.prmtop -r rec.prmtop -l lig.prmtop -n ":LIG" > ante_mmpbsa.log 2>&1 )
fi
cat > "$ANA/mmgbsa/mmgbsa.in" <<EOF
MM/GBSA $COMPUESTO-BRAF V600E (ML/MM $MODELO, $T_NS ns)
&general
  startframe=1, endframe=$N_FRAMES, interval=$INTERVALO, verbose=2, keep_files=0,
/
&gb
  igb=5, saltcon=0.150,
/
&decomp
  idecomp=1, dec_verbose=0,
/
EOF
MD_ABS=$(cd "$MD" && pwd)
DCDS=$(ls "$MD_ABS"/04_prod_*.dcd | sort | head -n "$SEGMENTOS" | tr '\n' ' ')
( cd "$ANA/mmgbsa" && en_md MMPBSA.py -O -i mmgbsa.in -o FINAL_MMGBSA.dat -do FINAL_DECOMP.dat \
      -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y $DCDS \
      -eo energias_por_frame.csv -deo decomp_por_frame.csv > mmpbsa.log 2>&1 )
en_md python 05_mmgbsa_figura.py --mmgbsa "$ANA/mmgbsa" --fig "$FIG" --dt-ps "$(calc "$INTERVALO * $GUARDAR_PS")"

echo "== 6) Comparación con la dinámica clásica ($CLASICO), primeros $T_NS ns"
if ls "$CLASICO"/04_prod_*.dcd > /dev/null 2>&1; then
    GB_C=""; [ -f "$GB_CLASICO/energias_por_frame.csv" ] && GB_C="--gb-clasico $GB_CLASICO"
    en_md python 06_comparar_clasico_ml.py --sistema ../sistema --clasico "$CLASICO" --ml "$MD" \
          $GB_C --gb-ml "$ANA/mmgbsa" --t-max "$T_NS" --out "$CMP" --titulo "$COMPUESTO"
else
    echo "   No hay corrida clásica en $CLASICO: correr antes 'bash correr_50ns.sh' para poder comparar."
fi

echo
echo "Listo. Trayectoria ML/MM: $MD   Análisis: $ANA   Figuras: $FIG   Comparación: $CMP"
