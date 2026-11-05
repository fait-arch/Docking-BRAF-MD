# -*- coding: utf-8 -*-
"""
Disposición de archivos de una ejecución.

El tiempo de simulación deja de estar codificado en nombres de carpeta sueltos
(``md_prueba``, ``md_50ns``, ``analisis_50ns``…) y pasa a ser una dimensión
explícita:

    simulaciones/<Compuesto>/
        sistema/                  <- compartido por todos los tiempos
        runs/
            10ns/{md,analisis,resultados,ml,estado}
            20ns/{md,analisis,resultados,ml,estado}
            50ns/{md,analisis,resultados,ml,estado}

Así el mismo compuesto se simula a 10, 20 y 50 ns sin que una ejecución pise a
la otra, y ambos pipelines leen la misma estructura.

**Compatibilidad**: si ya existen las carpetas antiguas de un tiempo dado
(``md_50ns``/``analisis_50ns``/``resultados_50ns``, o el trío de la prueba
rápida ``md_prueba``/``analisis``/``resultados_prueba``), el layout las reutiliza
en lugar de duplicar datos. Se puede forzar con ``modo='runs'`` o ``'legacy'``.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

#: Trío de carpetas heredado de la prueba rápida (100 ps).
_LEGACY_PRUEBA: Tuple[str, str, str] = ("md_prueba", "analisis", "resultados_prueba")

#: Tiempo (ns) al que corresponde la prueba rápida histórica.
NS_PRUEBA = 0.1


def etiqueta_tiempo(ns: float) -> str:
    """``50`` -> ``'50ns'``; ``0.1`` -> ``'0.1ns'``; ``20.0`` -> ``'20ns'``."""
    return f"{float(ns):g}ns"


def _legacy_candidato(etiqueta: str, ns: float) -> Tuple[str, str, str]:
    if abs(float(ns) - NS_PRUEBA) < 1e-9:
        return _LEGACY_PRUEBA
    return (f"md_{etiqueta}", f"analisis_{etiqueta}", f"resultados_{etiqueta}")


@dataclass
class RunLayout:
    """Rutas de una ejecución (compuesto + tiempo de simulación)."""

    base: Path          # directorio del compuesto
    ns: float           # tiempo de producción en ns
    etiqueta: str       # '50ns'
    md: Path
    analisis: Path
    resultados: Path
    ml: Path
    estado: Path
    legacy: bool        # True si reutiliza las carpetas antiguas

    # -- construcción -------------------------------------------------------
    @classmethod
    def resolver(cls, base: Path, ns: float, modo: str = "auto",
                 entradas: Optional[Path] = None) -> "RunLayout":
        """Construye el layout para ``(base, ns)``.

        ``modo``: ``'auto'`` reutiliza carpetas antiguas si existen, ``'runs'``
        fuerza la estructura nueva y ``'legacy'`` fuerza la antigua.
        """
        if modo not in ("auto", "runs", "legacy"):
            raise ValueError("modo debe ser 'auto', 'runs' o 'legacy'")
        base = Path(base)
        etiqueta = etiqueta_tiempo(ns)
        md_l, ana_l, res_l = _legacy_candidato(etiqueta, ns)
        raiz_nueva = base / "runs" / etiqueta

        usar_legacy = modo == "legacy"
        if modo == "auto":
            # Se reutiliza la carpeta antigua de MD salvo que ya exista una nueva
            # para ese mismo tiempo. Se mira 'runs/<t>/md' y no 'runs/<t>' entero,
            # porque 'runs/<t>/estado' se crea tambien en modo heredado y no debe
            # cambiar la resolucion en la siguiente ejecucion.
            usar_legacy = (base / md_l).exists() and not (raiz_nueva / "md").exists()

        if usar_legacy:
            # Se respetan las carpetas antiguas de MD y análisis (ahí están los
            # datos), pero lo nuevo —ML y estado— vive ya en runs/<tiempo>/,
            # para no sembrar más carpetas sueltas en la raíz del compuesto.
            return cls(
                base=base, ns=float(ns), etiqueta=etiqueta,
                md=base / md_l,
                analisis=base / ana_l,
                resultados=base / res_l,
                ml=raiz_nueva / "ml",
                estado=raiz_nueva / "estado",
                legacy=True,
                entradas=Path(entradas) if entradas else None,
            )
        return cls(
            base=base, ns=float(ns), etiqueta=etiqueta,
            md=raiz_nueva / "md",
            analisis=raiz_nueva / "analisis",
            resultados=raiz_nueva / "resultados",
            ml=raiz_nueva / "ml",
            estado=raiz_nueva / "estado",
            legacy=False,
            entradas=Path(entradas) if entradas else None,
        )

    # -- rutas compartidas --------------------------------------------------
    #: Carpeta de entradas del compuesto (scripts/, sistema/). Puede ser
    #: distinta de `base` cuando los resultados viven en otro árbol.
    entradas: Optional[Path] = None
    #: True cuando la corrida vive en una carpeta plana propia (10nsMacaridine/),
    #: donde la raíz de la corrida es `base` y no `base/runs/<tiempo>`.
    plano: bool = False

    @property
    def sistema(self) -> Path:
        """Topología AMBER: se prepara una vez y sirve para todos los tiempos.

        Vive con las entradas, no con los resultados: es la misma para 10, 20 o
        50 ns y es lo que garantiza que dos PC simulen el mismo sistema.
        """
        propio = self.base / "sistema"
        if self.entradas and not (propio / "sistema.prmtop").exists():
            return self.entradas / "sistema"
        return propio

    @property
    def raiz_run(self) -> Path:
        if self.plano:
            return self.base
        return self.base / "runs" / self.etiqueta

    @property
    def mmgbsa(self) -> Path:
        return self.analisis / "mmgbsa"

    @property
    def manifiesto(self) -> Path:
        return self.estado / "manifest.json"

    @property
    def log(self) -> Path:
        return self.estado / "ejecucion.log"

    @property
    def script(self) -> Path:
        return self.estado / "trabajo.sh"

    # -- utilidades ---------------------------------------------------------
    def crear(self) -> "RunLayout":
        """Crea los directorios necesarios."""
        for d in (self.md, self.analisis, self.resultados, self.ml, self.estado):
            d.mkdir(parents=True, exist_ok=True)
        return self

    def rel(self, ruta: Path) -> str:
        """Ruta relativa al directorio del compuesto (lo que esperan los scripts)."""
        try:
            return str(Path(ruta).resolve().relative_to(self.base.resolve()))
        except ValueError:
            return str(Path(ruta).resolve())

    def etapas_completadas(self) -> Dict[str, bool]:
        """Qué etapas de MD ya tienen su estado ``.xml`` guardado."""
        return {
            "min": (self.md / "01_min.xml").exists(),
            "cal": (self.md / "02_cal.xml").exists(),
            "eq": (self.md / "03_eq.xml").exists(),
        }

    def segmentos_hechos(self) -> int:
        return len(sorted(self.md.glob("04_prod_*.xml"))) if self.md.exists() else 0

    def trayectorias(self) -> list:
        return sorted(self.md.glob("04_prod_*.dcd")) if self.md.exists() else []

    def resumen_existente(self) -> Optional[Path]:
        for nombre in ("resumen.json", "reporte.json", "resultados.json"):
            p = self.analisis / nombre
            if p.exists():
                return p
        return None

    def __str__(self) -> str:  # pragma: no cover
        marca = " (carpetas heredadas)" if self.legacy else ""
        return f"{self.base.name} @ {self.etiqueta}{marca}"
