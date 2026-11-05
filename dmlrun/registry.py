# -*- coding: utf-8 -*-
"""
Registro de ejecuciones y manifiesto de reproducibilidad.

Cada envío deja dos rastros:

* una entrada en ``.dmlrun/registro.json`` (raíz del repo), que es lo que lee
  ``dmlrun status``;
* un ``manifest.json`` dentro de ``runs/<tiempo>/estado/`` con todo lo necesario
  para repetir la ejecución: commit de git, versiones de las bibliotecas que
  deciden el resultado, parámetros, semilla, servidor y plataforma.

El manifiesto es lo que permite comparar 10, 20 y 50 ns, o dos servidores
distintos, sabiendo exactamente qué cambió entre ellos. ``huella`` es un hash de
los parámetros científicos (sin el servidor ni las rutas): dos ejecuciones con la
misma huella son comparables aunque se hayan hecho en máquinas distintas.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import socket
import subprocess
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

DIR_REGISTRO = ".dmlrun"
ARCHIVO_REGISTRO = "registro.json"

#: Bibliotecas cuyas versiones afectan al resultado numérico.
_PAQUETES = ("openmm", "parmed", "MDAnalysis", "mdtraj", "numpy", "rdkit",
             "scikit-learn", "xgboost", "pandas")


def _ahora() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


def _git(raiz: Path, *args: str) -> Optional[str]:
    try:
        proc = subprocess.run(["git", "-C", str(raiz), *args],
                              capture_output=True, text=True, timeout=15)
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    return (proc.stdout or "").strip() or None


def versiones_paquetes() -> Dict[str, str]:
    """Versiones instaladas de las bibliotecas relevantes."""
    out: Dict[str, str] = {"python": platform.python_version()}
    try:
        from importlib.metadata import PackageNotFoundError, version
    except ImportError:  # pragma: no cover
        return out
    for nombre in _PAQUETES:
        try:
            out[nombre] = version(nombre)
        except PackageNotFoundError:
            continue
        except Exception:  # pragma: no cover
            continue
    return out


def herramientas_externas() -> Dict[str, Optional[str]]:
    """Rutas de los binarios externos que intervienen en el resultado."""
    import shutil
    return {h: shutil.which(h) for h in
            ("tleap", "antechamber", "parmchk2", "MMPBSA.py", "vina", "obabel",
             "pdb2pqr30", "gmx")}


@dataclass
class RunRecord:
    """Una ejecución registrada."""

    id: str
    compuesto: str
    pipeline: str
    tiempo_ns: float
    etiqueta: str
    servidor: str
    scheduler: str
    job_id: str
    estado: str
    creado: str
    directorio: str
    log: str
    script: str
    estado_dir: str
    etapas: List[str] = field(default_factory=list)
    huella: str = ""
    git_commit: Optional[str] = None
    notas: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, datos: Dict[str, Any]) -> "RunRecord":
        campos = set(cls.__dataclass_fields__)
        return cls(**{k: v for k, v in datos.items() if k in campos})


class Registro:
    """Lectura y escritura de ``.dmlrun/registro.json``."""

    def __init__(self, raiz: Path):
        self.raiz = Path(raiz)
        self.dir = self.raiz / DIR_REGISTRO
        self.archivo = self.dir / ARCHIVO_REGISTRO

    # -- lectura ------------------------------------------------------------
    def cargar(self) -> List[RunRecord]:
        if not self.archivo.exists():
            return []
        try:
            datos = json.loads(self.archivo.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return []
        registros = datos.get("ejecuciones", []) if isinstance(datos, dict) else datos
        out = []
        for d in registros:
            try:
                out.append(RunRecord.from_dict(d))
            except TypeError:
                continue
        return out

    def buscar(self, id_parcial: str) -> RunRecord:
        """Busca por id exacto o por prefijo."""
        registros = self.cargar()
        exactos = [r for r in registros if r.id == id_parcial]
        if exactos:
            return exactos[-1]
        parciales = [r for r in registros if r.id.startswith(id_parcial)]
        if len(parciales) == 1:
            return parciales[0]
        if not parciales:
            raise KeyError(f"no hay ninguna ejecución con id '{id_parcial}'")
        raise KeyError(
            f"'{id_parcial}' es ambiguo; coincide con: "
            + ", ".join(r.id for r in parciales[:6])
        )

    # -- escritura ----------------------------------------------------------
    def guardar(self, registros: List[RunRecord]) -> None:
        self.dir.mkdir(parents=True, exist_ok=True)
        payload = {
            "version": 1,
            "actualizado": _ahora(),
            "ejecuciones": [r.to_dict() for r in registros],
        }
        tmp = self.archivo.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False),
                       encoding="utf-8")
        os.replace(tmp, self.archivo)

    def añadir(self, registro: RunRecord) -> RunRecord:
        registros = self.cargar()
        registros.append(registro)
        self.guardar(registros)
        return registro

    def actualizar_estado(self, id_run: str, estado: str) -> None:
        registros = self.cargar()
        for r in registros:
            if r.id == id_run:
                r.estado = estado
        self.guardar(registros)


# ---------------------------------------------------------------------------
#  Manifiesto
# ---------------------------------------------------------------------------
def huella_parametros(datos: Dict[str, Any]) -> str:
    """Hash estable de los parámetros científicos (no del entorno)."""
    texto = json.dumps(datos, sort_keys=True, ensure_ascii=False, default=str)
    return hashlib.sha256(texto.encode("utf-8")).hexdigest()[:12]


def construir_manifiesto(*, raiz: Path, compuesto, layout, servidor,
                         pipeline, etapas: List[str],
                         parametros: Dict[str, Any]) -> Dict[str, Any]:
    """Manifiesto completo de una ejecución."""
    cientificos = {
        "compuesto": compuesto.clave,
        "pubchem": compuesto.pubchem,
        "smiles": compuesto.smiles,
        "carga_ligando": compuesto.carga_ligando,
        "receptor_residuos": compuesto.receptor_residuos,
        "ph": compuesto.ph,
        "pipeline": pipeline,
        "tiempo_ns": layout.ns,
        "parametros": parametros,
    }
    return {
        "version": 1,
        "creado": _ahora(),
        "huella": huella_parametros(cientificos),
        "cientifico": cientificos,
        "etapas": etapas,
        "entorno": {
            "servidor": servidor.nombre,
            "scheduler": servidor.scheduler,
            "plataforma_openmm": servidor.plataforma_openmm,
            "cpus": servidor.cpus,
            "gpus": servidor.gpus,
            "gres": servidor.gres,
            "hostname": socket.gethostname(),
            "sistema": platform.platform(),
            "conda_env": servidor.conda_env or os.environ.get("CONDA_DEFAULT_ENV"),
        },
        "codigo": {
            "git_commit": _git(raiz, "rev-parse", "HEAD"),
            "git_rama": _git(raiz, "rev-parse", "--abbrev-ref", "HEAD"),
            "git_sucio": bool(_git(raiz, "status", "--porcelain")),
            "raiz": str(raiz),
        },
        "paquetes": versiones_paquetes(),
        "herramientas": herramientas_externas(),
        "rutas": {
            "md": str(layout.md),
            "analisis": str(layout.analisis),
            "resultados": str(layout.resultados),
            "ml": str(layout.ml),
            "sistema": str(layout.sistema),
            "layout_heredado": layout.legacy,
        },
    }


def escribir_manifiesto(layout, manifiesto: Dict[str, Any]) -> Path:
    layout.estado.mkdir(parents=True, exist_ok=True)
    ruta = layout.manifiesto
    historico = layout.estado / "manifests"
    if ruta.exists():
        # no se pierde el manifiesto anterior: queda el histórico de la carpeta
        historico.mkdir(exist_ok=True)
        marca = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
        (historico / f"manifest-{marca}.json").write_text(
            ruta.read_text(encoding="utf-8"), encoding="utf-8")
    ruta.write_text(json.dumps(manifiesto, indent=2, ensure_ascii=False),
                    encoding="utf-8")
    return ruta
