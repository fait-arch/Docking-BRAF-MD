#!/usr/bin/env bash
# =============================================================================
#  Macaridine – BRAF V600E: simulación de 50 ns con los parámetros del artículo
#  (RSC Adv. 2022, 12, 30181): minimización 10 000 it. (k = 10 kcal/mol/Å²),
#  calentamiento 0→310 K en 500 ps, equilibrado 500 ps NPT, producción 50 ns NPT,
#  Langevin γ = 1 ps⁻¹, 1 atm, 2 fs, enlaces con H rígidos, corte 11 Å, PME.
#
#  Uso (desde esta carpeta, con el entorno activado):
#     bash correr_50ns.sh
#  Variables opcionales (ejemplos):
#     T_PROD_PS=100000 bash correr_50ns.sh      # 100 ns
#     REUSAR_MIN=1 bash correr_50ns.sh          # reutiliza md_prueba/01_min.* (mismo sistema)
#  Es reanudable: si se corta, al relanzarlo salta los segmentos de producción ya hechos.
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

T_CAL_PS=${T_CAL_PS:-500}          # calentamiento (ps)
T_EQ_PS=${T_EQ_PS:-500}            # equilibrado (ps)
T_PROD_PS=${T_PROD_PS:-50000}      # producción (ps) -> 50 ns
SEGMENTOS=${SEGMENTOS:-10}         # 10 segmentos de 5 ns (reanudables)
GUARDAR_PS=${GUARDAR_PS:-10}       # un fotograma cada 10 ps -> 5000 fotogramas
FRAMES_MMGBSA=${FRAMES_MMGBSA:-200}  # fotogramas usados en MM/GBSA
REUSAR_MIN=${REUSAR_MIN:-0}
BENCH_PASOS=${BENCH_PASOS:-2000}   # 0 = no medir velocidad
NOMBRE=${NOMBRE:-50ns}

MD=../md_${NOMBRE}
ANA=../analisis_${NOMBRE}
FIG=../resultados_${NOMBRE}
mkdir -p "$MD" "$ANA" "$FIG"

echo "== 0) Plataformas de OpenMM disponibles"
python - <<'PY'
import openmm as mm
nombres = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
print("   ", ", ".join(nombres))
if not any(n in nombres for n in ("CUDA", "HIP", "OpenCL")):
    print("    AVISO: no hay plataforma de GPU; en CPU 50 ns tardarían semanas.")
PY

echo "== 1) Velocidad esperada"
EST=../md_prueba/01_min.xml
if [ "$BENCH_PASOS" -gt 0 ] && [ -f "$EST" ]; then
    python benchmark.py --sistema ../sistema --estado "$EST" --pasos "$BENCH_PASOS" || true
fi

ETAPAS="min,cal,eq,prod"
if [ "$REUSAR_MIN" = "1" ] && [ -f ../md_prueba/01_min.xml ]; then
    cp ../md_prueba/01_min.xml ../md_prueba/01_min.pdb ../md_prueba/01_min.txt "$MD"/
    ETAPAS="cal,eq,prod"
fi
# si ya existen etapas terminadas, no repetirlas
[ -f "$MD/01_min.xml" ] && ETAPAS=${ETAPAS/min,/}
[ -f "$MD/02_cal.xml" ] && ETAPAS=${ETAPAS/cal,/}
[ -f "$MD/03_eq.xml" ]  && ETAPAS=${ETAPAS/eq,/}

echo "== 2) Dinámica molecular ($ETAPAS) -> $MD"
python 02_dinamica.py --sistema ../sistema --out "$MD" --etapas "$ETAPAS" \
       --t-calentamiento "$T_CAL_PS" --t-equilibrado "$T_EQ_PS" \
       --t-produccion "$T_PROD_PS" --segmentos "$SEGMENTOS" --guardar-cada "$GUARDAR_PS" \
       | tee -a "$MD/run.log"

echo "== 3) RMSD rápido"
python 03_rmsd.py --sistema ../sistema --md "$MD" --out "$FIG"

echo "== 4) Análisis completo y figuras"
python 04_analisis_reporte.py --base .. --md "$(basename "$MD")" --fig "$FIG" --out "$ANA"

echo "== 5) MM/GBSA ($FRAMES_MMGBSA fotogramas)"
mkdir -p "$ANA/mmgbsa"
N_FRAMES=$(python -c "print(int(round($T_PROD_PS / $GUARDAR_PS)))")
INTERVALO=$(( N_FRAMES / FRAMES_MMGBSA )); [ "$INTERVALO" -lt 1 ] && INTERVALO=1
if [ -f ../analisis/mmgbsa/com.prmtop ]; then
    cp ../analisis/mmgbsa/com.prmtop ../analisis/mmgbsa/rec.prmtop ../analisis/mmgbsa/lig.prmtop "$ANA/mmgbsa/"
else
    # topologías para MM/GBSA (radios mbondi2) generadas a partir del complejo seco
    ( cd "$ANA/mmgbsa" && rm -f com.prmtop rec.prmtop lig.prmtop \
      && printf 'changeRadii mbondi2\noutparm com.prmtop\nquit\n' \
         | parmed -n -p ../../sistema/complejo_seco.prmtop > parmed.log 2>&1 \
      && ante-MMPBSA.py -p com.prmtop -r rec.prmtop -l lig.prmtop -n ":LIG" > ante_mmpbsa.log 2>&1 )
fi
cat > "$ANA/mmgbsa/mmgbsa.in" <<EOF
MM/GBSA Macaridine-BRAF V600E ($NOMBRE)
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
DCDS=$(ls "$MD_ABS"/04_prod_*.dcd | sort | tr '\n' ' ')
( cd "$ANA/mmgbsa" && MMPBSA.py -O -i mmgbsa.in -o FINAL_MMGBSA.dat -do FINAL_DECOMP.dat \
      -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y $DCDS \
      -eo energias_por_frame.csv -deo decomp_por_frame.csv > mmpbsa.log 2>&1 )
python 05_mmgbsa_figura.py --mmgbsa "$ANA/mmgbsa" --fig "$FIG" \
       --dt-ps "$(python -c "print($INTERVALO * $GUARDAR_PS)")"

echo
echo "Listo. Trayectoria: $MD   Análisis: $ANA   Figuras: $FIG"
