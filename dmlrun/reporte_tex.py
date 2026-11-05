# -*- coding: utf-8 -*-
"""
Generación de las secciones LaTeX del reporte, dentro del propio documento.

El documento (``paper/simulacion.tex``) está escrito a mano y no debe tocarse:
este módulo **solo** reemplaza lo que hay entre sus propias marcas

    % <<<dmlrun:nombre>>>
    ... contenido generado ...
    % <<<dmlrun:fin:nombre>>>

y deja intacto todo lo demás. Si las marcas no existen, se insertan una vez
justo antes de ``\\appendix``. Así el reporte se puede regenerar tantas veces
como haga falta sin perder una sola línea de la redacción del autor.

Se generan dos bloques:

``datos_50ns``
    Tablas con las cifras de cada corrida de 50 ns, leídas de ``resumen.json`` y
    de ``FINAL_MMGBSA.dat``. Al generarse del dato y no copiarse a mano, no se
    desincronizan del análisis.

``comparacion_md_ml``
    La comparación que pide el trabajo: el mismo conjunto de compuestos
    ordenado por la vía tradicional (docking + dinámica de 50 ns + MM/GBSA) y
    por la vía de aprendizaje automático (docking + rescoring), con la
    discusión de si coinciden.
"""
from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

MARCA_INI = "% <<<dmlrun:{nombre}>>>"
MARCA_FIN = "% <<<dmlrun:fin:{nombre}>>>"

#: Dónde se insertan los bloques nuevos si aún no existen sus marcas: la línea
#: que contiene únicamente \appendix. Se exige la línea entera porque buscar la
#: cadena suelta tambien casa dentro de \appendixname, \appendixtitle y
#: cualquier otra macro que empiece igual, y partiria el preámbulo por la mitad.
ANCLA_INSERCION = re.compile(r"(?m)^[ \t]*\\appendix[ \t]*$")
ANCLA_CIERRE = re.compile(r"(?m)^[ \t]*\\end\{document\}")


# ---------------------------------------------------------------------------
#  Utilidades de LaTeX
# ---------------------------------------------------------------------------
_ESCAPES = {
    "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
    "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
}


def tex_escape(texto) -> str:
    """Texto plano -> LaTeX seguro."""
    s = str(texto)
    s = s.replace("\\", r"\textbackslash{}")
    for a, b in _ESCAPES.items():
        s = s.replace(a, b)
    return s


def num(valor, dec: int = 2) -> str:
    """Número para siunitx, con los decimales que no anulen la magnitud."""
    if valor is None:
        return "--"
    try:
        v = float(valor)
    except (TypeError, ValueError):
        return tex_escape(valor)
    for d in (dec, 3, 4):
        texto = f"{v:.{d}f}"
        if float(texto) != 0.0 or v == 0.0:
            return texto
    return f"{v:.2e}"


def media_sd(valor, dec: int = 2) -> str:
    """``[media, sd, ...]`` -> ``\\num{1.95 +- 0.28}`` en notación del documento."""
    if not isinstance(valor, (list, tuple)) or not valor:
        return "--"
    if len(valor) >= 2:
        return f"${num(valor[0], dec)} \\pm {num(valor[1], dec)}$"
    return f"${num(valor[0], dec)}$"


def atomos_pesados(formula: Optional[str]) -> int:
    """Cuenta los átomos distintos de hidrógeno de una fórmula molecular."""
    if not formula:
        return 0
    pares = re.findall(r"([A-Z][a-z]?)(\d*)", str(formula))
    return sum(int(c or 1) for el, c in pares if el and el != "H")


# ---------------------------------------------------------------------------
#  Lectura de los resultados de una corrida
# ---------------------------------------------------------------------------
@dataclass
class Corrida:
    """Las cifras de un sistema a un tiempo dado."""

    clave: str
    nombre: str
    etiqueta: str
    ns: float
    resumen: Dict
    mmgbsa: Optional[Tuple[float, float]] = None       # (media, DE)
    mmgbsa_componentes: Optional[Dict[str, Tuple[float, float]]] = None
    ml: Optional[Dict] = None                          # score del modelo, si lo hay

    # -- derivados ----------------------------------------------------------
    @property
    def vina(self) -> Optional[float]:
        v = (self.resumen.get("vina") or {}).get("mejor")
        return float(v) if v is not None else None

    @property
    def formula(self) -> Optional[str]:
        return (self.resumen.get("ligando") or {}).get("formula")

    @property
    def n_pesados(self) -> int:
        return atomos_pesados(self.formula)

    @property
    def eficiencia(self) -> Optional[float]:
        """ΔG MM/GBSA por átomo pesado."""
        if self.mmgbsa and self.n_pesados:
            return self.mmgbsa[0] / self.n_pesados
        return None

    @property
    def eficiencia_vina(self) -> Optional[float]:
        if self.vina is not None and self.n_pesados:
            return self.vina / self.n_pesados
        return None

    @property
    def rmsd_ligando(self) -> Optional[Sequence[float]]:
        return self.resumen.get("rmsd_lig_prod")

    @property
    def pose_estable(self) -> Optional[bool]:
        """¿Sobrevivió la pose del docking a la dinámica?

        Criterio: RMSD medio del ligando por debajo de 3 Å respecto a la pose
        acoplada. Por encima, el compuesto se ha reacomodado y la pose que
        puntuó el docking ya no es la que se simula.
        """
        r = self.rmsd_ligando
        if not isinstance(r, (list, tuple)) or not r:
            return None
        return float(r[0]) < 3.0


def _parsear_mmgbsa(dat: Path) -> Tuple[Optional[Tuple[float, float]], Dict]:
    """ΔG total y componentes de FINAL_MMGBSA.dat."""
    try:
        texto = Path(dat).read_text(errors="replace")
    except OSError:
        return None, {}
    etiquetas = {
        "VDWAALS": "van der Waals",
        "EEL": "Electrostática",
        "EGB": "Solvatación polar",
        "ESURF": "Solvatación no polar",
    }
    comp: Dict[str, Tuple[float, float]] = {}
    bloque = texto.split("DELTA Energy Terms")[-1]
    for linea in bloque.splitlines():
        for clave, nombre in etiquetas.items():
            if linea.strip().startswith(clave):
                nums = re.findall(r"-?\d+\.\d+", linea)
                if len(nums) >= 2:
                    comp[nombre] = (float(nums[0]), float(nums[1]))
                break
    m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", texto)
    total = (float(m.group(1)), float(m.group(2))) if m else None
    return total, comp


def _leer_ml(ml_dir: Path) -> Optional[Dict]:
    """Mejor puntuación del rescoring, si el pipeline de ML ya se ejecutó."""
    csvf = Path(ml_dir) / "rescored.csv"
    if not csvf.exists():
        return None
    try:
        with open(csvf, newline="", encoding="utf-8") as fh:
            filas = list(csv.DictReader(fh))
    except OSError:
        return None
    if not filas:
        return None
    columna = next((c for c in ("final_score", "ml_score", "score")
                    if c in filas[0]), None)
    if not columna:
        return None
    valores = []
    for f in filas:
        try:
            valores.append(float(f[columna]))
        except (TypeError, ValueError):
            continue
    if not valores:
        return None
    return {"columna": columna, "mejor": max(valores),
            "media": sum(valores) / len(valores), "n": len(valores)}


def cargar_corrida(compuesto, layout) -> Optional[Corrida]:
    """Lee los resultados de un sistema a un tiempo; ``None`` si no hay análisis."""
    resumen_path = Path(layout.analisis) / "resumen.json"
    if not resumen_path.exists():
        return None
    try:
        resumen = json.loads(resumen_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    total, comp = _parsear_mmgbsa(Path(layout.mmgbsa) / "FINAL_MMGBSA.dat")
    return Corrida(
        clave=compuesto.clave, nombre=compuesto.nombre,
        etiqueta=layout.etiqueta, ns=float(layout.ns),
        resumen=resumen, mmgbsa=total, mmgbsa_componentes=comp,
        ml=_leer_ml(layout.ml),
    )


# ---------------------------------------------------------------------------
#  Bloque 1: tablas de datos de las corridas
# ---------------------------------------------------------------------------
def bloque_datos(corridas: List[Corrida]) -> str:
    if not corridas:
        return "% (sin corridas con análisis todavía)\n"

    filas_est = []
    for c in corridas:
        r = c.resumen
        filas_est.append(
            f"{tex_escape(c.nombre)} & "
            f"{media_sd(r.get('rmsd_CA_prod'))} & "
            f"{media_sd(r.get('rmsd_lig_prod'))} & "
            f"{media_sd(r.get('rg_prod'))} & "
            f"{media_sd(r.get('ligando_enterrado_pct'), 1)} \\\\"
        )

    filas_ene = []
    for c in corridas:
        mm = (f"${num(c.mmgbsa[0], 1)} \\pm {num(c.mmgbsa[1], 1)}$"
              if c.mmgbsa else "--")
        ef = f"${num(c.eficiencia, 2)}$" if c.eficiencia is not None else "--"
        filas_ene.append(
            f"{tex_escape(c.nombre)} & "
            f"{num(c.vina, 2) if c.vina is not None else '--'} & "
            f"{mm} & {c.n_pesados or '--'} & {ef} \\\\"
        )

    # componentes de MM/GBSA, una columna por sistema
    nombres = [tex_escape(c.nombre) for c in corridas]
    componentes = ["van der Waals", "Electrostática",
                   "Solvatación polar", "Solvatación no polar"]
    filas_comp = []
    for nombre_c in componentes:
        celdas = []
        for c in corridas:
            val = (c.mmgbsa_componentes or {}).get(nombre_c)
            celdas.append(f"${num(val[0], 1)} \\pm {num(val[1], 1)}$" if val else "--")
        filas_comp.append(f"{nombre_c} & " + " & ".join(celdas) + " \\\\")
    total_celdas = []
    for c in corridas:
        total_celdas.append(f"${num(c.mmgbsa[0], 1)} \\pm {num(c.mmgbsa[1], 1)}$"
                            if c.mmgbsa else "--")
    filas_comp.append("\\midrule")
    filas_comp.append("$\\Delta G$ total & " + " & ".join(total_celdas) + " \\\\")

    cols_comp = "l" + "c" * len(corridas)
    nl = "\n"
    return f"""% Generado por `dmlrun reporte --formato tex`; se regenera con el mismo comando.
% No editar a mano: lo que haya entre las marcas se reemplaza.
\\section{{Resumen cuantitativo de las corridas de \\SI{{50}}{{ns}}}}
\\label{{sec:datos50ns}}

Las cifras de esta sección se leen directamente de los archivos de análisis de
cada corrida (\\texttt{{analisis\\_50ns/resumen.json}} y
\\texttt{{FINAL\\_MMGBSA.dat}}), de modo que no pueden quedar desfasadas respecto
a los datos.

\\begin{{table}}[H]
\\centering\\small
\\begin{{tabular}}{{@{{}}lcccc@{{}}}}
\\toprule
Sistema & RMSD C$\\alpha$ (\\AA) & RMSD ligando (\\AA) & $R_g$ (\\AA) & Enterrado (\\%) \\\\
\\midrule
{nl.join(filas_est)}
\\bottomrule
\\end{{tabular}}
\\caption{{Estabilidad estructural durante la producción (media $\\pm$ desviación
estándar). El RMSD del ligando se mide respecto a la pose de partida: valores
por encima de \\SI{{3}}{{\\angstrom}} indican que el compuesto se reacomodó y que
la pose puntuada por el \\textit{{docking}} no es la que se simuló.}}
\\label{{tab:estabilidad50}}
\\end{{table}}

\\begin{{table}}[H]
\\centering\\small
\\begin{{tabular}}{{@{{}}lcccc@{{}}}}
\\toprule
Sistema & Vina (kcal/mol) & $\\Delta G$ MM/GBSA (kcal/mol) & Átomos pesados & $\\Delta G$/átomo \\\\
\\midrule
{nl.join(filas_ene)}
\\bottomrule
\\end{{tabular}}
\\caption{{Energías de unión por los dos métodos. La última columna es la
eficiencia del ligando: normaliza el $\\Delta G$ por el tamaño de la molécula y
evita confundir una afinidad mayor con un mejor candidato cuando los compuestos
tienen distinto número de átomos. El MM/GBSA no incluye el término entrópico
$-T\\Delta S$, así que sobrestima la afinidad en valor absoluto y solo debe
usarse para ordenar sistemas calculados con el mismo protocolo.}}
\\label{{tab:energias50}}
\\end{{table}}

\\begin{{table}}[H]
\\centering\\small
\\begin{{tabular}}{{@{{}}{cols_comp}@{{}}}}
\\toprule
Término (kcal/mol) & {" & ".join(nombres)} \\\\
\\midrule
{nl.join(filas_comp)}
\\bottomrule
\\end{{tabular}}
\\caption{{Descomposición del MM/GBSA. Separa la contribución hidrofóbica
(van der Waals) de la electrostática, que es la que cambia al protonar el
ligando.}}
\\label{{tab:mmgbsacomp}}
\\end{{table}}
"""


# ---------------------------------------------------------------------------
#  Bloque 2: comparación entre la vía tradicional y la de aprendizaje automático
# ---------------------------------------------------------------------------
def _orden(corridas: List[Corrida], clave, descendente=False) -> List[str]:
    """Orden de los sistemas según una magnitud; ``None`` va al final."""
    con_valor = [c for c in corridas if clave(c) is not None]
    sin_valor = [c for c in corridas if clave(c) is None]
    ordenados = sorted(con_valor, key=clave, reverse=descendente)
    return [c.nombre for c in ordenados + sin_valor]


def bloque_comparacion(corridas: List[Corrida]) -> str:
    if not corridas:
        return "% (sin corridas con análisis todavía)\n"

    hay_ml = any(c.ml for c in corridas)

    # --- tabla de rankings -------------------------------------------------
    # vía rápida: el docking ordena de más negativo a menos
    orden_vina = _orden(corridas, lambda c: c.vina)
    # vía tradicional: MM/GBSA tras la dinámica
    orden_mm = _orden(corridas, lambda c: c.mmgbsa[0] if c.mmgbsa else None)
    # vía de aprendizaje automático: mayor puntuación = mejor
    orden_ml = (_orden(corridas, lambda c: c.ml["mejor"] if c.ml else None,
                       descendente=True) if hay_ml else [])

    filas = []
    for i, c in enumerate(corridas, start=1):
        puesto = lambda orden: (str(orden.index(c.nombre) + 1)
                                if c.nombre in orden else "--")
        ml_celda = (f"${num(c.ml['mejor'], 3)}$" if c.ml else "--")
        estable = {True: "sí", False: "no", None: "--"}[c.pose_estable]
        filas.append(
            f"{tex_escape(c.nombre)} & "
            f"{num(c.vina, 2) if c.vina is not None else '--'} ({puesto(orden_vina)}) & "
            + (f"{ml_celda} ({puesto(orden_ml)}) & " if hay_ml else "")
            + (f"${num(c.mmgbsa[0], 1)}$ ({puesto(orden_mm)})" if c.mmgbsa else "-- (--)")
            + f" & {estable} \\\\"
        )

    cols = "lcc c" if not hay_ml else "lccc c"
    cab_ml = "Rescoring ML (puesto) & " if hay_ml else ""
    nl = "\n"

    # --- lectura del resultado --------------------------------------------
    coinciden = (orden_vina == orden_mm) if (orden_vina and orden_mm) else None
    inestables = [c for c in corridas if c.pose_estable is False]

    # El orden coincidente solo significa algo si las poses comparadas son las
    # que el sistema mantiene: si el mejor clasificado se desplaza durante la
    # produccion, el acuerdo entre metodos es superficial, porque cada uno esta
    # describiendo una geometria distinta.
    primero_inestable = None
    if orden_mm:
        primero_inestable = next((c for c in corridas
                                  if c.nombre == orden_mm[0]), None)
    top_inestable = bool(primero_inestable and primero_inestable.pose_estable is False)

    if coinciden is True and not top_inestable:
        lectura = ("Los dos métodos ordenan los compuestos igual, y las poses se "
                   "mantienen durante la producción. Cuando ocurren las dos cosas "
                   "a la vez, el cribado rápido habría bastado para priorizar, y "
                   "la dinámica aporta el mecanismo —qué contactos sostienen la "
                   "unión— más que el orden.")
    elif coinciden is True and top_inestable:
        lectura = (f"Los dos métodos ordenan igual los compuestos, pero "
                   f"\\textbf{{el acuerdo es engañoso}}: "
                   f"{tex_escape(primero_inestable.nombre)}, que ambos sitúan "
                   "primero, no conserva la pose acoplada durante la producción. "
                   "El \\textit{docking} puntuó una geometría y la dinámica "
                   "acabó describiendo otra, de modo que la coincidencia de "
                   "puestos no indica que el cribado rápido fuese suficiente. "
                   "El MM/GBSA de un ligando que se desplaza tampoco es fiable "
                   "como afinidad: promedia conformaciones que ya no son la que "
                   "se quería evaluar.")
    elif coinciden is False:
        mejor_rapido = orden_vina[0] if orden_vina else "--"
        mejor_lento = orden_mm[0] if orden_mm else "--"
        lectura = (f"Los dos métodos \\textbf{{no coinciden}}: el "
                   f"\\textit{{docking}} sitúa primero a {tex_escape(mejor_rapido)} "
                   f"y el MM/GBSA tras la dinámica a {tex_escape(mejor_lento)}. "
                   "Es el caso que justifica el coste de la dinámica: la "
                   "discrepancia no es ruido del método, sino que el "
                   "\\textit{docking} puntúa una pose fija y la dinámica "
                   "comprueba si esa pose sobrevive.")
    else:
        lectura = ("Falta todavía alguno de los dos resultados para poder "
                   "comparar los órdenes.")

    if inestables:
        nombres_i = ", ".join(tex_escape(c.nombre) for c in inestables)
        detalle_estabilidad = (
            f"La última columna es la que lo detecta. En {nombres_i} el ligando se "
            f"desplaza más de \\SI{{3}}{{\\angstrom}} respecto a la pose acoplada, "
            "de modo que la energía que había calculado el \\textit{docking} "
            "corresponde a una geometría que el sistema abandona. Ninguna "
            "función de puntuación sobre una pose estática —ni clásica ni "
            "aprendida— puede detectar eso, porque la información no está en la "
            "pose sino en su evolución temporal.")
    else:
        detalle_estabilidad = (
            "En todos los sistemas el ligando se mantiene cerca de la pose "
            "acoplada durante la producción, así que las puntuaciones sobre la "
            "pose estática describen la misma geometría que la dinámica.")

    if hay_ml:
        bloque_ml = (
            "El rescoring con aprendizaje automático se calculó sobre las mismas "
            "poses, con descriptores de contacto por capas (OnionNet) y huellas "
            "de interacción, y un modelo entrenado aparte. Comparte con el "
            "\\textit{docking} la limitación de partir de una estructura fija, "
            "pero sustituye una función de puntuación fijada a mano por una "
            "ajustada a datos de afinidad.")
        pendiente = ""
    else:
        bloque_ml = (
            "La columna de aprendizaje automático está vacía porque todavía no "
            "hay un modelo entrenado: el rescoring necesita un conjunto de "
            "compuestos con afinidad medida frente a BRAF (por ejemplo, de "
            "ChEMBL) para ajustar el modelo antes de poder puntuar estos.")
        pendiente = r"""
\begin{quote}\small
\textbf{Qué falta para completar esta comparación.} Entrenar el modelo con
afinidades conocidas y volver a generar esta sección:
\begin{verbatim}
dml-train -X descriptores.csv -y afinidades.csv -m xgb -o modelos/rescorer.pkl
dmlrun run --pipeline ml --compuesto todos --tiempo 50 \
       --fuente docking --modelo modelos/rescorer.pkl --servidor hpc-cpu
dmlrun reporte --formato tex
\end{verbatim}
Mientras tanto, la tabla compara lo que sí está calculado: la función de
puntuación clásica del \textit{docking} frente al MM/GBSA sobre la dinámica.
\end{quote}
"""

    return f"""% Generado por `dmlrun reporte --formato tex`; se regenera con el mismo comando.
% No editar a mano: lo que haya entre las marcas se reemplaza.
\\section{{Comparación: dinámica tradicional frente a aprendizaje automático}}
\\label{{sec:mdvsml}}

Las dos vías responden a la misma pregunta —¿cuál de estos compuestos se une
mejor al sitio de ATP de BRAF V600E?— con costes muy distintos. La vía
tradicional acopla, simula \\SI{{50}}{{ns}} y estima la energía con MM/GBSA sobre
la trayectoria. La vía de aprendizaje automático acopla y vuelve a puntuar las
poses con un modelo entrenado, sin simular. Esta sección compara lo que cada
una concluye sobre \\emph{{el mismo}} conjunto de compuestos.

\\begin{{table}}[H]
\\centering\\small
\\begin{{tabular}}{{@{{}}{cols}@{{}}}}
\\toprule
Sistema & Vina (puesto) & {cab_ml}$\\Delta G$ MM/GBSA (puesto) & Pose estable \\\\
\\midrule
{nl.join(filas)}
\\bottomrule
\\end{{tabular}}
\\caption{{El mismo conjunto de compuestos ordenado por cada método. Entre
paréntesis, el puesto en el \\textit{{ranking}} de ese método. ``Pose estable''
indica si el RMSD medio del ligando durante la producción se mantuvo por debajo
de \\SI{{3}}{{\\angstrom}} respecto a la pose de partida.}}
\\label{{tab:mdvsml}}
\\end{{table}}

\\subsection{{Qué dice la comparación}}

{lectura}

{detalle_estabilidad}

{bloque_ml}
{pendiente}
\\subsection{{Qué método usar}}

Los dos métodos no compiten: responden preguntas distintas y el orden natural es
usarlos en cascada.

\\begin{{itemize}}
\\item El \\textit{{docking}} y el rescoring con aprendizaje automático son
  métodos de \\textbf{{cribado}}: puntúan una pose fija, cuestan segundos o
  minutos por compuesto y sirven para reducir miles de candidatos a decenas. Su
  punto ciego es que no saben si la pose que puntúan es la que el compuesto
  adopta realmente.
\\item La dinámica molecular con MM/GBSA es un método de \\textbf{{verificación}}:
  cuesta horas de GPU por compuesto, así que no se aplica a miles, pero es la
  que revela si la pose se mantiene, qué contactos persisten y cuáles se pierden.
\\item Por eso la comparación de la tabla~\\ref{{tab:mdvsml}} no busca un ganador,
  sino medir \\textbf{{cuánto se pierde al saltarse la verificación}}. Cuando los
  órdenes coinciden, el cribado basta para priorizar; cuando no, la dinámica
  está detectando algo que ninguna puntuación sobre estructura fija podía ver.
\\end{{itemize}}

Un modelo de aprendizaje automático entrenado sobre descriptores de la
trayectoria, y no solo de la pose acoplada, sería el paso siguiente: conservaría
parte de la información que aporta la dinámica a un coste mucho menor, aunque
seguiría necesitando dinámicas previas para generar sus datos de entrenamiento.
"""


# ---------------------------------------------------------------------------
#  Inserción en el documento
# ---------------------------------------------------------------------------
def insertar_bloque(texto: str, nombre: str, contenido: str) -> Tuple[str, str]:
    """Reemplaza el bloque entre marcas; lo crea si no existe.

    Devuelve ``(texto_nuevo, accion)`` donde acción es ``'reemplazado'`` o
    ``'insertado'``. Nunca toca nada fuera de las marcas.
    """
    ini = MARCA_INI.format(nombre=nombre)
    fin = MARCA_FIN.format(nombre=nombre)
    bloque = f"{ini}\n{contenido.rstrip()}\n{fin}"

    if ini in texto and fin in texto:
        inicio = texto.index(ini)
        final = texto.index(fin) + len(fin)
        return texto[:inicio] + bloque + texto[final:], "reemplazado"

    if ini in texto or fin in texto:
        raise ValueError(
            f"el documento tiene solo una de las dos marcas de '{nombre}'; "
            "corrígelo a mano antes de volver a generar")

    m = ANCLA_INSERCION.search(texto)
    if m:
        pos = m.start()
        # delante del ancla, separado por líneas en blanco
        return texto[:pos] + bloque + "\n\n\n" + texto[pos:], "insertado"

    m = ANCLA_CIERRE.search(texto)
    if m:
        pos = m.start()
        return texto[:pos] + bloque + "\n\n" + texto[pos:], "insertado"
    return texto.rstrip() + "\n\n" + bloque + "\n", "insertado"


def escribir_en_documento(documento: Path, bloques: Dict[str, str],
                          copia_seguridad: bool = True) -> Dict[str, str]:
    """Aplica todos los bloques al documento y devuelve qué se hizo con cada uno."""
    documento = Path(documento)
    texto = documento.read_text(encoding="utf-8")
    original = texto
    acciones: Dict[str, str] = {}
    for nombre, contenido in bloques.items():
        texto, accion = insertar_bloque(texto, nombre, contenido)
        acciones[nombre] = accion
    if texto == original:
        return acciones
    if copia_seguridad:
        marca = datetime.now().strftime("%Y%m%d-%H%M%S")
        respaldo = documento.with_name(f"{documento.stem}.{marca}.bak")
        respaldo.write_text(original, encoding="utf-8")
        acciones["_respaldo"] = str(respaldo)
    documento.write_text(texto, encoding="utf-8")
    return acciones
