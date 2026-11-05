# -*- coding: utf-8 -*-
"""
Recolección de resultados y generación del informe.

Las figuras y las tablas de una corrida quedan repartidas entre la carpeta de
resultados, la de análisis y (en las corridas antiguas) ``Img/<Compuesto>``.
Este módulo las reúne en una sola carpeta por corrida y escribe un
``informe.html`` autocontenido: las imágenes van incrustadas, así que el archivo
se puede abrir sin conexión, copiar a un pendrive o enviar por correo y se ve
igual en cualquier ordenador.

El informe está escrito para que lo lea alguien que no ha lanzado la simulación:
cada figura lleva qué muestra y **qué hay que mirar en ella**, y los números van
en tablas con sus unidades. Al lado quedan los PDF de cada figura (para
publicación) y los CSV (para rehacer cualquier gráfica).
"""
from __future__ import annotations

import base64
import csv
import html
import json
import re
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
#  Catálogo de figuras: nombre de archivo -> título, qué es y qué mirar
# ---------------------------------------------------------------------------
FIGURAS: List[Dict[str, str]] = [
    {
        "archivo": "fig_ligando_2d",
        "titulo": "Estructura del ligando",
        "que_es": "Estructura 2D del compuesto tal como se parametrizó, con la "
                  "protonación usada en la simulación.",
        "que_mirar": "Comprobar que es el compuesto correcto y que el estado de "
                     "protonación (cargas, grupos OH/NH) es el que se esperaba a "
                     "pH 7.4. Si la amina aparece neutra y en el organismo estaría "
                     "protonada, los contactos electrostáticos saldrán subestimados.",
    },
    {
        "archivo": "fig_vina",
        "titulo": "Energías de acoplamiento (AutoDock Vina)",
        "que_es": "Distribución de las energías de las poses generadas por el "
                  "docking, antes de la dinámica.",
        "que_mirar": "La mejor pose (más negativa) es la que entra en la dinámica. "
                     "Un grupo de poses con energías parecidas indica que el sitio "
                     "admite varias orientaciones; conviene entonces desconfiar de "
                     "una sola pose y mirar qué hace el ligando durante la MD.",
    },
    {
        "archivo": "fig_complejo_general",
        "titulo": "Vista general del complejo",
        "que_es": "Render del complejo proteína-ligando completo.",
        "que_mirar": "Situar el ligando dentro de la quinasa y confirmar que está "
                     "en el sitio de ATP y no en una cavidad superficial.",
    },
    {
        "archivo": "fig_bolsillo",
        "titulo": "Bolsillo de unión",
        "que_es": "Detalle del sitio de ATP con los residuos que rodean al ligando.",
        "que_mirar": "Qué residuos envuelven al ligando y si hay huecos sin ocupar. "
                     "Es la figura que explica, en términos estructurales, los "
                     "contactos que luego aparecen cuantificados.",
    },
    {
        "archivo": "fig_termodinamica",
        "titulo": "Control del sistema: temperatura, densidad, energía y volumen",
        "que_es": "Las cuatro magnitudes de control a lo largo de las tres fases "
                  "(calentamiento, equilibrado y producción).",
        "que_mirar": "**Es la figura de control de calidad: si falla, lo demás no "
                     "es interpretable.** La temperatura debe estabilizarse en "
                     "310 K y la densidad en torno a 1.01 g/cm³, ambas planas "
                     "durante la producción. Una deriva sostenida de la energía "
                     "potencial o del volumen indica que el sistema no estaba "
                     "equilibrado cuando empezó la producción.",
    },
    {
        "archivo": "fig_rmsd",
        "titulo": "Estabilidad: RMSD de la proteína y del ligando",
        "que_es": "Desviación respecto a la estructura inicial, por separado para "
                  "los carbonos α de la proteína y para el ligando.",
        "que_mirar": "La curva de Cα debe alcanzar una meseta (típicamente por "
                     "debajo de 2-3 Å): significa que la proteína ya no cambia de "
                     "conformación y la parte estable de la trayectoria es la que "
                     "se analiza. El RMSD del ligando es la clave farmacológica: "
                     "si se estabiliza, el compuesto se queda donde se acopló; si "
                     "sube de forma sostenida, se está reacomodando o saliendo del "
                     "sitio, y entonces la pose del docking no se sostiene.",
    },
    {
        "archivo": "fig_rg",
        "titulo": "Radio de giro",
        "que_es": "Compactación global de la proteína durante la trayectoria.",
        "que_mirar": "Debe mantenerse prácticamente constante. Un aumento "
                     "progresivo indicaría que el dominio se está abriendo o "
                     "desplegando, lo que invalidaría el análisis del sitio activo.",
    },
    {
        "archivo": "fig_distancias",
        "titulo": "Distancias a los residuos clave del sitio de ATP",
        "que_es": "Distancia del ligando a Lys483, Cys532 (bisagra) y Asp594 "
                  "(motivo DFG) a lo largo del tiempo.",
        "que_mirar": "Son los tres anclajes que definen un inhibidor de BRAF. "
                     "Distancias estables y cortas (< 4-5 Å) indican un modo de "
                     "unión mantenido; si una se abre y no vuelve, el ligando ha "
                     "perdido ese contacto. Lo interesante es cuál pierde y cuál "
                     "conserva, porque distingue inhibidores de tipo I y II.",
    },
    {
        "archivo": "fig_rmsf",
        "titulo": "Flexibilidad por residuo (RMSF)",
        "que_es": "Fluctuación de cada residuo, con las regiones funcionales de la "
                  "quinasa señaladas (lazo P, hélice αC, bisagra, lazo catalítico, "
                  "DFG y lazo de activación).",
        "que_mirar": "Los picos marcan las zonas móviles. Lo esperable es que los "
                     "extremos y el lazo de activación sean los más flexibles y "
                     "que el lazo catalítico y la bisagra estén rígidos. Una "
                     "rigidez anómala en el lazo de activación junto al ligando "
                     "sugiere que el compuesto lo está bloqueando, que es "
                     "precisamente el mecanismo buscado.",
    },
    {
        "archivo": "fig_contactos",
        "titulo": "Contactos ligando-proteína y su persistencia",
        "que_es": "Residuos que contactan con el ligando y fracción del tiempo de "
                  "producción en que se mantiene cada contacto.",
        "que_mirar": "**Es la figura que define el farmacóforo real**, no el del "
                     "docking. Los residuos con ocupación cercana a 1 están en "
                     "contacto casi todo el tiempo y son los que habría que "
                     "conservar al modificar el compuesto o mutar para validar "
                     "experimentalmente. Los de ocupación baja son contactos "
                     "ocasionales y no deberían guiar el diseño.",
    },
    {
        "archivo": "fig_mmgbsa",
        "titulo": "Energía libre de unión (MM/GBSA) y su descomposición",
        "que_es": "ΔG de unión estimado por MM/GBSA a lo largo de la trayectoria y "
                  "contribución de cada residuo.",
        "que_mirar": "La media con su desviación es la estimación de afinidad; "
                     "la descomposición dice **qué residuos la producen**. Dos "
                     "advertencias al interpretarla: MM/GBSA no incluye el término "
                     "entrópico −TΔS, así que los valores son más negativos que la "
                     "afinidad real y sólo sirven para comparar compuestos entre sí "
                     "con el mismo protocolo; y un ΔG más favorable por mayor "
                     "tamaño del ligando no significa mejor candidato, para eso "
                     "está la eficiencia por átomo pesado de la tabla.",
    },
]

#: Tablas de resultados: clave en resumen.json -> cómo presentarla
ESTADISTICOS = [
    ("rmsd_CA_prod", "RMSD de Cα (producción)", "Å",
     "Estabilidad de la proteína. Meseta baja = trayectoria aprovechable."),
    ("rmsd_lig_prod", "RMSD del ligando (producción)", "Å",
     "Si es bajo y estable, el ligando se mantiene en la pose acoplada."),
    ("rg_prod", "Radio de giro", "Å",
     "Compactación global; debe ser constante."),
    ("d_K483_prod", "Distancia a Lys483", "Å",
     "Anclaje catalítico del sitio de ATP."),
    ("d_C532_prod", "Distancia a Cys532 (bisagra)", "Å",
     "Contacto de bisagra, característico de los inhibidores de BRAF."),
    ("d_D594_prod", "Distancia a Asp594 (DFG)", "Å",
     "Contacto con el motivo DFG; informa del tipo de inhibición."),
    ("T_prod", "Temperatura", "K", "Control: debe ser 310 K."),
    ("rho_prod", "Densidad", "g/cm³", "Control: debe ser ≈ 1.01."),
]

LIGANDO = [
    ("formula", "Fórmula", ""),
    ("MW", "Masa molecular", "g/mol"),
    ("logP", "logP", ""),
    ("HBD", "Donadores de puente de H", ""),
    ("HBA", "Aceptores de puente de H", ""),
    ("TPSA", "Superficie polar (TPSA)", "Å²"),
    ("rotables", "Enlaces rotables", ""),
]


# ---------------------------------------------------------------------------
#  Recolección
# ---------------------------------------------------------------------------
@dataclass
class Recoleccion:
    """Lo que se encontró de una corrida."""

    compuesto: object
    layout: object
    figuras: Dict[str, Dict[str, Path]] = field(default_factory=dict)  # base -> {png, pdf}
    datos: List[Path] = field(default_factory=list)
    resumen: Optional[Dict] = None
    mmgbsa: Optional[Dict] = None
    manifiesto: Optional[Dict] = None

    @property
    def vacia(self) -> bool:
        return not self.figuras and not self.resumen


def _leer_json(ruta: Path) -> Optional[Dict]:
    try:
        return json.loads(Path(ruta).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, ValueError):
        return None


def _parsear_mmgbsa(dat: Path) -> Optional[Dict]:
    """Extrae ΔG total y sus componentes de FINAL_MMGBSA.dat."""
    try:
        texto = Path(dat).read_text(errors="replace")
    except OSError:
        return None
    out: Dict[str, Tuple[float, float]] = {}
    etiquetas = {
        "DELTA TOTAL": "ΔG total",
        "VDWAALS": "van der Waals",
        "EEL": "Electrostática",
        "EGB": "Solvatación polar",
        "ESURF": "Solvatación no polar",
    }
    # el bloque de diferencias es el último del archivo
    bloque = texto.split("DELTA Energy Terms")[-1]
    for linea in bloque.splitlines():
        for clave, nombre in etiquetas.items():
            if linea.strip().startswith(clave):
                nums = re.findall(r"-?\d+\.\d+", linea)
                if len(nums) >= 2:
                    out[nombre] = (float(nums[0]), float(nums[1]))
                break
    m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", texto)
    if m and "ΔG total" not in out:
        out["ΔG total"] = (float(m.group(1)), float(m.group(2)))
    return out or None


def recolectar(compuesto, layout, raiz_img: Optional[Path] = None) -> Recoleccion:
    """Busca figuras y datos de una corrida en todos los sitios posibles."""
    rec = Recoleccion(compuesto=compuesto, layout=layout)

    # Carpetas donde pueden estar las figuras, de la más específica a la más
    # antigua. La carpeta global Img/<Compuesto> sólo se mira en las corridas
    # heredadas, porque es donde iban las figuras antes de que el tiempo de
    # simulación fuese una carpeta propia; mirarla siempre haría que una corrida
    # nueva y vacía apareciese con las figuras de otra.
    carpetas = [layout.resultados, layout.analisis]
    # Img/<Compuesto> no lleva el tiempo en el nombre: guarda las figuras de la
    # ultima corrida que se analizo. Solo se mira si esta corrida tiene su propio
    # analisis (resumen.json); si no, las figuras de ahi serian de otra corrida y
    # apareceria un informe con graficas que no le corresponden.
    tiene_analisis = (Path(layout.analisis) / "resumen.json").exists()
    if raiz_img and getattr(layout, "legacy", False) and tiene_analisis:
        carpetas.append(Path(raiz_img) / compuesto.directorio.name)

    for carpeta in carpetas:
        if not Path(carpeta).is_dir():
            continue
        for img in sorted(Path(carpeta).glob("*.png")) + sorted(Path(carpeta).glob("*.pdf")):
            base = img.stem
            entrada = rec.figuras.setdefault(base, {})
            # la primera carpeta que la tiene gana (la más específica)
            entrada.setdefault(img.suffix.lstrip("."), img)

    for carpeta in (layout.analisis, layout.resultados):
        if Path(carpeta).is_dir():
            rec.datos += sorted(Path(carpeta).glob("*.csv"))

    rec.resumen = _leer_json(layout.analisis / "resumen.json")
    dat = layout.mmgbsa / "FINAL_MMGBSA.dat"
    if dat.exists():
        rec.mmgbsa = _parsear_mmgbsa(dat)
    rec.manifiesto = _leer_json(layout.manifiesto)
    return rec


# ---------------------------------------------------------------------------
#  HTML
# ---------------------------------------------------------------------------
CSS = """
:root {
  --sup: #fcfcfb; --tarjeta: #ffffff; --tinta: #16150f; --tinta2: #5c5a52;
  --linea: #e4e3df; --azul: #2a78d6; --naranja: #eb6834; --aviso: #fff8e6;
  --aviso-linea: #e8c97a;
}
* { box-sizing: border-box; }
body {
  margin: 0; padding: 0 16px 72px; background: var(--sup); color: var(--tinta);
  font: 16px/1.65 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
        "Helvetica Neue", Arial, sans-serif;
  -webkit-font-smoothing: antialiased;
}
.envoltura { max-width: 940px; margin: 0 auto; }
header { padding: 56px 0 28px; border-bottom: 2px solid var(--tinta); margin-bottom: 36px; }
h1 { font-size: 1.95rem; line-height: 1.2; margin: 0 0 6px; letter-spacing: -0.015em; }
.sub { color: var(--tinta2); font-size: 1.02rem; margin: 0; }
h2 {
  font-size: 1.3rem; margin: 48px 0 6px; padding-top: 22px;
  border-top: 1px solid var(--linea); letter-spacing: -0.01em;
}
h2:first-of-type { border-top: 0; padding-top: 0; }
h3 { font-size: 1.06rem; margin: 30px 0 4px; }
p { margin: 0 0 14px; }
.meta { display: flex; flex-wrap: wrap; gap: 8px 28px; margin: 18px 0 0;
        font-size: 0.86rem; color: var(--tinta2); }
.meta b { color: var(--tinta); font-weight: 600; }
figure { margin: 18px 0 8px; }
figure img {
  width: 100%; height: auto; display: block; border: 1px solid var(--linea);
  border-radius: 8px; background: #fff;
}
figcaption { font-size: 0.9rem; color: var(--tinta2); margin-top: 10px; }
.mirar {
  background: var(--tarjeta); border: 1px solid var(--linea);
  border-left: 3px solid var(--azul); border-radius: 0 8px 8px 0;
  padding: 13px 16px; margin: 12px 0 0; font-size: 0.92rem;
}
.mirar b.et { display: block; font-size: 0.74rem; letter-spacing: 0.08em;
              text-transform: uppercase; color: var(--azul); margin-bottom: 5px; }
table { width: 100%; border-collapse: collapse; margin: 16px 0 8px; font-size: 0.92rem; }
th, td { text-align: left; padding: 9px 12px; border-bottom: 1px solid var(--linea); }
th { font-weight: 600; font-size: 0.78rem; letter-spacing: 0.05em;
     text-transform: uppercase; color: var(--tinta2); }
td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
tr:last-child td { border-bottom: 0; }
.nota { color: var(--tinta2); font-size: 0.88rem; }
.aviso {
  background: var(--aviso); border: 1px solid var(--aviso-linea);
  border-radius: 8px; padding: 14px 16px; margin: 20px 0; font-size: 0.93rem;
}
.aviso b { display: block; margin-bottom: 4px; }
.archivos { font-size: 0.9rem; }
.archivos code { background: #f3f2ee; padding: 1px 5px; border-radius: 4px;
                 font-size: 0.86em; }
footer { margin-top: 56px; padding-top: 20px; border-top: 1px solid var(--linea);
         font-size: 0.84rem; color: var(--tinta2); }
@media print {
  body { padding: 0; font-size: 11pt; }
  h2 { page-break-after: avoid; } figure { page-break-inside: avoid; }
  .mirar { border-left-color: #888; }
}
"""


def _b64(ruta: Path) -> Optional[str]:
    try:
        return base64.b64encode(Path(ruta).read_bytes()).decode("ascii")
    except OSError:
        return None


def _fmt(v, dec: int = 2) -> str:
    """Numero legible, con decimales suficientes para no perder la magnitud.

    Una desviacion como 0.0025 g/cm3 no puede imprimirse como "0": si el
    redondeo normal la anula, se aumentan los decimales hasta que se vea.
    """
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return html.escape(str(v))
    if isinstance(v, int):
        return str(v)
    for d in (dec, 3, 4, 5):
        texto = f"{v:.{d}f}"
        if float(texto) != 0.0 or v == 0.0:
            return texto.rstrip("0").rstrip(".") if "." in texto else texto
    return f"{v:.2e}"


def _stat(valor) -> str:
    """``[media, sd, min, max]`` -> ``'1.39 ± 0.08  (1.21 – 1.56)'``."""
    if not isinstance(valor, (list, tuple)) or not valor:
        return _fmt(valor)
    if len(valor) >= 4:
        return (f"{_fmt(valor[0])} ± {_fmt(valor[1])}"
                f" <span class='nota'>({_fmt(valor[2])} – {_fmt(valor[3])})</span>")
    if len(valor) >= 2:
        return f"{_fmt(valor[0])} ± {_fmt(valor[1])}"
    return _fmt(valor[0])


def _tabla(filas: List[Tuple[str, str, str]], cabecera: Tuple[str, str, str]) -> str:
    if not filas:
        return ""
    th = "".join(f"<th>{html.escape(c)}</th>" for c in cabecera)
    cuerpo = ""
    for etiqueta, valor, nota in filas:
        cuerpo += (f"<tr><td>{etiqueta}</td><td class='num'>{valor}</td>"
                   f"<td class='nota'>{nota}</td></tr>")
    return f"<table><thead><tr>{th}</tr></thead><tbody>{cuerpo}</tbody></table>"


def _seccion_resumen(rec: Recoleccion) -> str:
    r = rec.resumen or {}
    if not r:
        return ""
    out = ["<h2>Resultados numéricos</h2>"]

    filas = []
    for clave, etiqueta, unidad, nota in ESTADISTICOS:
        if clave not in r:
            continue
        valor = _stat(r[clave])
        filas.append((f"{html.escape(etiqueta)}"
                      + (f" <span class='nota'>[{html.escape(unidad)}]</span>" if unidad else ""),
                      valor, nota))
    if filas:
        out.append("<h3>Estabilidad y geometría</h3>")
        out.append(_tabla(filas, ("Magnitud", "Media ± DE (mín – máx)", "Qué indica")))

    # contactos
    contactos = r.get("contactos_50") or []
    if contactos:
        filas = [(f"{html.escape(str(c[1]))}{c[0]}", f"{float(c[2]) * 100:.0f} %",
                  "contacto mantenido" if float(c[2]) > 0.9 else "contacto parcial")
                 for c in contactos[:12]]
        out.append("<h3>Residuos en contacto con el ligando</h3>")
        out.append("<p class='nota'>Fracción del tiempo de producción en que se "
                   "mantiene el contacto. Los que superan el 90 % son los "
                   "candidatos a mutagénesis dirigida.</p>")
        out.append(_tabla(filas, ("Residuo", "Ocupación", "Lectura")))

    # puentes de hidrógeno
    ph = r.get("puentes_H")
    frac = r.get("frames_con_puenteH")
    if ph is not None:
        if ph:
            # cada entrada es [donador, aceptor, ocupación]; el último valor es
            # la fracción de tiempo y los anteriores identifican los átomos
            filas = []
            for p in ph[:10]:
                partes = list(p) if isinstance(p, (list, tuple)) else [p]
                try:
                    ocupacion = float(partes[-1])
                except (TypeError, ValueError):
                    continue
                par = " → ".join(html.escape(str(x)) for x in partes[:-1])
                lectura = ("contacto permanente" if ocupacion > 0.9
                           else "contacto intermitente")
                filas.append((par or "—", f"{ocupacion * 100:.0f} %", lectura))
            out.append("<h3>Puentes de hidrógeno</h3>")
            out.append("<p class='nota'>Un puente de hidrógeno con ocupación alta "
                       "es el contacto más específico que puede tener el compuesto: "
                       "orienta la molécula y suele explicar la selectividad frente "
                       "a otras quinasas.</p>")
            out.append(_tabla(filas, ("Donador → aceptor", "Ocupación", "Lectura")))
        else:
            out.append("<h3>Puentes de hidrógeno</h3>")
            out.append("<p>No se detectó ningún puente de hidrógeno estable entre "
                       "el ligando y la proteína"
                       + (f" (fotogramas con al menos uno: {float(frac) * 100:.0f} %)."
                          if frac is not None else ".")
                       + " La unión sería entonces mayoritariamente hidrofóbica, lo "
                       "que suele traducirse en menor selectividad.</p>")

    # RMSF por región
    reg = r.get("rmsf_regiones") or {}
    if reg:
        filas = [(html.escape(k), f"{_fmt(v)} Å", "") for k, v in reg.items()]
        out.append("<h3>Flexibilidad por región funcional</h3>")
        out.append(_tabla(filas, ("Región", "RMSF medio", "")))

    # enterramiento y docking
    filas = []
    ent = r.get("ligando_enterrado_pct")
    if isinstance(ent, (list, tuple)) and ent:
        filas.append(("Ligando enterrado", f"{_fmt(ent[0])} %",
                      "Fracción de la superficie del ligando oculta por la proteína; "
                      "valores altos indican que permanece dentro del bolsillo."))
    vina = r.get("vina") or {}
    if vina.get("mejor") is not None:
        filas.append(("Mejor pose de docking (Vina)", f"{_fmt(vina['mejor'])} kcal/mol",
                      f"De {vina.get('n', '?')} poses generadas."))
    if filas:
        out.append("<h3>Acoplamiento y enterramiento</h3>")
        out.append(_tabla(filas, ("Magnitud", "Valor", "Qué indica")))

    # ligando
    lig = r.get("ligando") or {}
    if lig:
        filas = [(html.escape(etiqueta), f"{_fmt(lig[clave])}"
                  + (f" {html.escape(unidad)}" if unidad else ""), "")
                 for clave, etiqueta, unidad in LIGANDO if clave in lig]
        out.append("<h3>Propiedades del compuesto</h3>")
        out.append(_tabla(filas, ("Propiedad", "Valor", "")))

    return "\n".join(out)


def _seccion_mmgbsa(rec: Recoleccion) -> str:
    if not rec.mmgbsa:
        return ""
    r = rec.resumen or {}
    filas = []
    for nombre, (media, sd) in rec.mmgbsa.items():
        filas.append((html.escape(nombre), f"{media:.1f} ± {sd:.1f}", ""))
    out = ["<h3>Energía libre de unión (MM/GBSA)</h3>",
           _tabla(filas, ("Término", "kcal/mol", ""))]
    total = rec.mmgbsa.get("ΔG total")
    pesados, n = None, 0
    lig = r.get("ligando") or {}
    formula = str(lig.get("formula") or "")
    if total and formula:
        # átomos pesados = todos los de la fórmula menos los hidrógenos
        pares = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
        n = sum(int(c or 1) for el, c in pares if el and el != "H")
        if n:
            pesados = total[0] / n
    if pesados is not None:
        out.append(f"<p class='nota'>Eficiencia por átomo pesado: "
                   f"<b>{pesados:.2f} kcal/mol por átomo</b> "
                   f"({n} átomos pesados). Es la cifra a comparar entre compuestos "
                   f"de tamaño distinto: un ΔG más negativo sólo por ser una "
                   f"molécula más grande no indica un mejor candidato.</p>")
    out.append("<div class='aviso'><b>Cómo leer este número</b>"
               "El cálculo no incluye el término entrópico −TΔS, así que sobrestima "
               "la afinidad en términos absolutos. Sirve para ordenar compuestos "
               "calculados con el mismo protocolo, no para predecir una constante de "
               "disociación.</div>")
    return "\n".join(out)


def _aviso_duracion(rec: Recoleccion) -> str:
    ns = float(rec.layout.ns)
    if ns >= 20:
        return ""
    if ns < 1:
        return ("<div class='aviso'><b>Esta corrida es una prueba técnica, no un "
                f"resultado</b>La producción duró {ns * 1000:.0f} ps. Sirve para "
                "comprobar que el sistema está bien montado y para medir la "
                "velocidad, pero en ese tiempo una proteína no explora sus "
                "conformaciones: ni el RMSD, ni los contactos, ni el MM/GBSA deben "
                "interpretarse como propiedades del complejo. Para conclusiones hay "
                "que repetirlo a 20-50 ns o más.</div>")
    return ("<div class='aviso'><b>Duración corta</b>"
            f"La producción duró {ns:g} ns. Es suficiente para ver si el ligando se "
            "mantiene en el sitio, pero corta para afirmar que el complejo está "
            "equilibrado o para comparar afinidades con confianza. Conviene "
            "contrastarla con una corrida más larga del mismo compuesto.</div>")


def _seccion_figuras(rec: Recoleccion, carpeta_figuras: Path) -> str:
    out: List[str] = ["<h2>Figuras</h2>"]
    puestas = set()
    for spec in FIGURAS:
        base = spec["archivo"]
        # coincidencia exacta o por prefijo (fig_bolsillo_A, fig_bolsillo_B…)
        candidatos = sorted(k for k in rec.figuras
                            if k == base or k.startswith(base + "_"))
        for clave in candidatos:
            if clave in puestas:
                continue
            fuentes = rec.figuras[clave]
            png = fuentes.get("png")
            pdf = fuentes.get("pdf")
            if not png and not pdf:
                continue
            puestas.add(clave)
            titulo = spec["titulo"]
            if clave != base:
                titulo += f" — {clave.replace(base + '_', '')}"
            out.append(f"<h3>{html.escape(titulo)}</h3>")
            if png:
                datos = _b64(png)
                if datos:
                    out.append(f'<figure><img alt="{html.escape(titulo)}" '
                               f'src="data:image/png;base64,{datos}">')
                    out.append(f"<figcaption>{spec['que_es']}</figcaption></figure>")
            formatos = []
            for ext, archivo in sorted(fuentes.items()):
                destino = carpeta_figuras / f"{clave}.{ext}"
                formatos.append(f"<code>figuras/{destino.name}</code>")
            out.append(f"<div class='mirar'><b class='et'>Qué mirar</b>"
                       f"{spec['que_mirar']}</div>")
            if formatos:
                out.append(f"<p class='archivos nota'>Archivo: "
                           f"{', '.join(formatos)}</p>")

    # figuras encontradas que no están en el catálogo
    resto = [k for k in sorted(rec.figuras) if k not in puestas]
    # rmsd.png (del script rapido) repite fig_rmsd, que es la version completa
    if "fig_rmsd" in rec.figuras:
        resto = [k for k in resto if k != "rmsd"]
    if resto:
        out.append("<h3>Otras figuras</h3>")
        for clave in resto:
            png = rec.figuras[clave].get("png")
            if not png:
                continue
            datos = _b64(png)
            if datos:
                out.append(f'<figure><img alt="{html.escape(clave)}" '
                           f'src="data:image/png;base64,{datos}">'
                           f"<figcaption>{html.escape(clave)}</figcaption></figure>")
    return "\n".join(out)


def generar_html(rec: Recoleccion, carpeta: Path) -> str:
    """Informe completo de una corrida, autocontenido."""
    C, L = rec.compuesto, rec.layout
    r = rec.resumen or {}
    man = rec.manifiesto or {}
    dur = r.get("duraciones_ps") or {}
    prod_ps = dur.get("produccion")

    meta = [
        ("Compuesto", f"{C.nombre}" + (f" (PubChem {C.pubchem})" if C.pubchem else "")),
        ("Diana", "BRAF V600E, dominio quinasa " + str(C.receptor_residuos)),
        ("Producción", f"{prod_ps / 1000:g} ns" if prod_ps else f"{L.ns:g} ns"),
        ("Protocolo", "310 K, 1 atm, Langevin γ = 1 ps⁻¹, paso 2 fs, PME, corte 11 Å"),
    ]
    if dur.get("calentamiento"):
        meta.append(("Equilibración",
                     f"{dur['calentamiento'] / 1000:g} ns calentamiento + "
                     f"{dur.get('equilibrado', 0) / 1000:g} ns NPT"))
    ent = (man.get("entorno") or {})
    if ent.get("servidor"):
        meta.append(("Calculado en", f"{ent['servidor']} ({ent.get('plataforma_openmm', '-')})"))
    if man.get("huella"):
        meta.append(("Huella de parámetros", man["huella"]))
    cod = man.get("codigo") or {}
    if cod.get("git_commit"):
        meta.append(("Versión del código", str(cod["git_commit"])[:10]))

    meta_html = "".join(f"<span><b>{html.escape(k)}:</b> {html.escape(str(v))}</span>"
                        for k, v in meta)

    partes = [
        f"<header><h1>{html.escape(C.nombre)} – BRAF V600E</h1>",
        f"<p class='sub'>Dinámica molecular y análisis de la unión · "
        f"{html.escape(L.etiqueta)}</p>",
        f"<div class='meta'>{meta_html}</div></header>",
        _aviso_duracion(rec),
        _seccion_figuras(rec, carpeta / "figuras"),
        _seccion_resumen(rec),
        _seccion_mmgbsa(rec),
    ]

    if rec.datos:
        nombres = ", ".join(f"<code>datos/{p.name}</code>" for p in rec.datos[:12])
        partes.append(
            "<h2>Datos para rehacer las gráficas</h2>"
            f"<p class='archivos'>Las series completas están en {nombres}. "
            "Cada figura se puede regenerar o reestilizar desde su CSV, y los PDF "
            "de <code>figuras/</code> son vectoriales, válidos para publicación.</p>")

    partes.append(
        "<footer>Informe generado con <code>dmlrun reporte</code> el "
        f"{datetime.now().strftime('%d/%m/%Y a las %H:%M')}. "
        "Las imágenes van incrustadas: este archivo se puede abrir sin conexión, "
        "copiar o enviar por correo tal cual.</footer>")

    cuerpo = "\n".join(p for p in partes if p)
    return f"""<!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(C.nombre)} – BRAF V600E · {html.escape(L.etiqueta)}</title>
<style>{CSS}</style>
</head>
<body><div class="envoltura">
{cuerpo}
</div></body>
</html>
"""


# ---------------------------------------------------------------------------
#  Escritura de la carpeta del informe
# ---------------------------------------------------------------------------
def escribir_reporte(rec: Recoleccion, destino: Path) -> Dict[str, object]:
    """Crea ``destino`` con informe.html, figuras/ y datos/."""
    destino = Path(destino)
    figs = destino / "figuras"
    dats = destino / "datos"
    for d in (destino, figs, dats):
        d.mkdir(parents=True, exist_ok=True)

    copiadas = 0
    for clave, fuentes in rec.figuras.items():
        for ext, archivo in fuentes.items():
            try:
                shutil.copy2(archivo, figs / f"{clave}.{ext}")
                copiadas += 1
            except OSError:
                continue
    for csvf in rec.datos:
        try:
            shutil.copy2(csvf, dats / csvf.name)
        except OSError:
            continue
    for extra in (rec.layout.analisis / "resumen.json",
                  rec.layout.mmgbsa / "FINAL_MMGBSA.dat",
                  rec.layout.mmgbsa / "FINAL_DECOMP.dat",
                  rec.layout.manifiesto):
        if Path(extra).exists():
            try:
                shutil.copy2(extra, dats / Path(extra).name)
            except OSError:
                continue

    informe = destino / "informe.html"
    informe.write_text(generar_html(rec, destino), encoding="utf-8")
    return {"informe": informe, "figuras": copiadas, "datos": len(rec.datos),
            "carpeta": destino}


# ---------------------------------------------------------------------------
#  Índice y tabla comparativa
# ---------------------------------------------------------------------------
def _valor_comparable(rec: Recoleccion, clave: str) -> str:
    r = rec.resumen or {}
    if clave == "mmgbsa":
        if rec.mmgbsa and "ΔG total" in rec.mmgbsa:
            m, s = rec.mmgbsa["ΔG total"]
            return f"{m:.1f} ± {s:.1f}"
        return "—"
    v = r.get(clave)
    if isinstance(v, (list, tuple)) and len(v) >= 2:
        return f"{_fmt(v[0])} ± {_fmt(v[1])}"
    if v is None:
        return "—"
    return _fmt(v)


def escribir_indice(recolecciones: List[Tuple[Recoleccion, Path]],
                    destino: Path) -> Path:
    """Índice con enlaces a cada informe y una tabla comparativa."""
    filas = ""
    for rec, carpeta in recolecciones:
        r = rec.resumen or {}
        ent = ((rec.manifiesto or {}).get("entorno") or {})
        rel = f"{carpeta.name}/informe.html"
        enterrado = r.get("ligando_enterrado_pct")
        filas += (
            "<tr>"
            f"<td><a href='{html.escape(rel)}'>{html.escape(rec.compuesto.nombre)}</a></td>"
            f"<td>{html.escape(rec.layout.etiqueta)}</td>"
            f"<td class='num'>{_valor_comparable(rec, 'rmsd_CA_prod')}</td>"
            f"<td class='num'>{_valor_comparable(rec, 'rmsd_lig_prod')}</td>"
            f"<td class='num'>{_valor_comparable(rec, 'mmgbsa')}</td>"
            f"<td class='num'>"
            f"{(_fmt(enterrado[0]) + ' %') if isinstance(enterrado, (list, tuple)) and enterrado else '—'}"
            f"</td>"
            f"<td>{html.escape(str(ent.get('servidor', '—')))}</td>"
            "</tr>")

    # Un indice sin contexto invita a leer como resultados lo que son pruebas
    # tecnicas: se avisa de las corridas demasiado cortas antes de la tabla.
    cortas = [f"{r.compuesto.nombre} ({r.layout.etiqueta})"
              for r, _ in recolecciones if float(r.layout.ns) < 20]
    aviso = ""
    if cortas:
        aviso = ("<div class='aviso'><b>Hay corridas demasiado cortas en esta "
                 "tabla</b>" + ", ".join(html.escape(c) for c in cortas) +
                 ". En ese tiempo la proteina no explora sus conformaciones, asi "
                 "que sus cifras sirven para comprobar que el montaje funciona, "
                 "no como propiedades del complejo. Cada informe lo indica "
                 "tambien en su cabecera.</div>")

    html_doc = f"""<!doctype html>
<html lang="es">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Resultados – BRAF V600E</title>
<style>{CSS}</style>
</head>
<body><div class="envoltura">
<header>
  <h1>Resultados de dinámica molecular</h1>
  <p class="sub">BRAF V600E · un informe por compuesto y tiempo de simulación</p>
</header>

{aviso}
<h2>Comparación</h2>
<p class='nota'>Haz clic en el compuesto para abrir su informe completo con las
figuras y su interpretación. El RMSD del ligando es el indicador de si se
mantiene en el sitio de unión; el MM/GBSA sólo es comparable entre filas
calculadas con el mismo protocolo y sin el término entrópico.</p>
<table>
<thead><tr>
  <th>Compuesto</th><th>Tiempo</th><th>RMSD Cα (Å)</th><th>RMSD ligando (Å)</th>
  <th>ΔG MM/GBSA (kcal/mol)</th><th>Enterrado</th><th>Servidor</th>
</tr></thead>
<tbody>{filas}</tbody>
</table>

<footer>Generado con <code>dmlrun reporte</code> el
{datetime.now().strftime('%d/%m/%Y a las %H:%M')}. Carpeta autocontenida: se
puede copiar o comprimir entera y abrir en cualquier ordenador sin conexión.</footer>
</div></body>
</html>
"""
    destino = Path(destino)
    destino.parent.mkdir(parents=True, exist_ok=True)
    destino.write_text(html_doc, encoding="utf-8")
    return destino
