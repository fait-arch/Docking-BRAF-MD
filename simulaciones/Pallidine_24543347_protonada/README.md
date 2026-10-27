# Pallidine protonada (PubChem 24543347, NH⁺) + BRAF V600E

Variante de `../Pallidine_24543347` con la **amina terciaria protonada** (catión, carga +1), la forma
probable a pH 7.4 (p*K*a de las morfinandienonas ≈ 8). Todo lo demás es idéntico: receptor, protocolo,
parámetros y scripts.

## Cómo se preparó

- **Misma pose de partida que la forma neutra** (RMSD 0.000 Å en átomos pesados): se tomó
  `../Pallidine_24543347/sistema/ligando_pose.sdf` (mejor pose de Vina, −8.07 kcal/mol) y se añadió el H⁺ al
  N (d N–H = 1.03 Å), en la dirección de su par libre. Así la única diferencia entre las dos simulaciones es
  la protonación. El N queda como nuevo estereocentro (*R*); los carbonos conservan su configuración (*S*,*S*).
- `entrada/Pallidine24543347_protonada_ligand.pdb`: ese ligando protonado (46 átomos).
- Sistema:
  ```bash
  python 01_preparar_sistema.py --receptor ../entrada/Pallidine24543347_receptor.pdb \
         --ligando ../entrada/Pallidine24543347_protonada_ligand.pdb --out ../sistema \
         --sin-docking --carga-ligando 1
  ```
  GAFF + AM1-BCC con carga +1 (q(N) = −0.68, q(H–N) = +0.45). Receptor idéntico al de la forma neutra.
  Carga del soluto +9 → SLTCAP: 36 Na⁺ + 45 Cl⁻ (la forma neutra: 37 + 45). 49 097 átomos, octaedro
  a = 89.51 Å.
- En la pose inicial, el N⁺ está a 4.2 Å del carboxilato de Asp594 (motivo DFG): puede formar un puente salino.

## Correr

```bash
cd scripts
python 02_dinamica.py --sistema ../sistema --out ../md_prueba --etapas min
python benchmark.py --sistema ../sistema --estado ../md_prueba/01_min.xml --pasos 2000
REUSAR_MIN=1 nohup bash correr_50ns.sh > corrida_50ns.log 2>&1 &
```

Resultados en `md_50ns/`, `analisis_50ns/` y `resultados_50ns/` (ignorados por git). Comparar con
`../Pallidine_24543347/` (misma corrida con el ligando neutro).
