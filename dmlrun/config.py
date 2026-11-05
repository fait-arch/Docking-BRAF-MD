# -*- coding: utf-8 -*-
"""
Carga de la configuración declarativa del proyecto.

Dos archivos en ``config/``:

* ``servidores.yaml``  -> perfiles de entorno de ejecución (``PerfilServidor``).
* ``compuestos.yaml``  -> compuestos y tiempos predefinidos (``Compuesto``).

Nada de esto se codifica en los pipelines: la CLI resuelve el perfil a partir
del nombre que se pasa con ``--servidor`` y lo inyecta en el trabajo generado.
"""
from __future__ import annotations

import os
import shutil
import socket
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

ARCHIVO_SERVIDORES = "config/servidores.yaml"
ARCHIVO_COMPUESTOS = "config/compuestos.yaml"

#: Marcadores que identifican la raíz del repositorio.
_MARCADORES = (ARCHIVO_SERVIDORES, "dockml", "setup.py")


class ConfigError(RuntimeError):
    """Error de configuración legible para el usuario."""


def _cargar_yaml(ruta: Path) -> Dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover
        raise ConfigError(
            "Falta PyYAML. Instálalo en el entorno:\n"
            "    conda install -c conda-forge pyyaml   # o: pip install pyyaml"
        ) from exc
    if not ruta.exists():
        raise ConfigError(f"No existe el archivo de configuración: {ruta}")
    with open(ruta, "r", encoding="utf-8") as fh:
        datos = yaml.safe_load(fh) or {}
    if not isinstance(datos, dict):
        raise ConfigError(f"{ruta} debe contener un mapa YAML en la raíz")
    return datos


def encontrar_raiz(inicio: Optional[os.PathLike] = None) -> Path:
    """Devuelve la raíz del repositorio dockingML.

    Orden de búsqueda: ``$DMLRUN_RAIZ``, el argumento ``inicio``, el directorio
    actual y sus padres, y por último la ubicación de este propio paquete.
    """
    env = os.environ.get("DMLRUN_RAIZ")
    if env:
        raiz = Path(env).expanduser().resolve()
        if (raiz / ARCHIVO_SERVIDORES).exists():
            return raiz
        raise ConfigError(
            f"$DMLRUN_RAIZ={raiz} no contiene {ARCHIVO_SERVIDORES}"
        )

    candidatos: List[Path] = []
    if inicio is not None:
        candidatos.append(Path(inicio).expanduser().resolve())
    candidatos.append(Path.cwd().resolve())
    candidatos.append(Path(__file__).resolve().parent.parent)

    for base in candidatos:
        for d in [base, *base.parents]:
            if (d / ARCHIVO_SERVIDORES).exists():
                return d
    # Último recurso: un directorio que tenga el paquete dockml.
    for base in candidatos:
        for d in [base, *base.parents]:
            if (d / "dockml").is_dir() and (d / "setup.py").exists():
                return d
    raise ConfigError(
        "No encuentro la raíz del proyecto. Ejecuta la CLI dentro del repositorio "
        "o exporta DMLRUN_RAIZ=/ruta/a/dockingML"
    )


# ---------------------------------------------------------------------------
#  Perfiles de servidor
# ---------------------------------------------------------------------------
@dataclass
class PerfilServidor:
    """Entorno de ejecución resuelto: scheduler, recursos y preparación."""

    nombre: str
    descripcion: str = ""
    scheduler: str = "local"
    cpus: int = 4
    memoria_gb: int = 16
    horas: int = 8
    gpus: int = 0
    gres: Optional[str] = None
    plataforma_openmm: str = "auto"
    conda_env: Optional[str] = None
    #: Como se activa el entorno: 'micromamba' (este proyecto) o 'conda'.
    gestor_entorno: str = "micromamba"
    modulos: List[str] = field(default_factory=list)
    pre_comandos: List[str] = field(default_factory=list)
    post_comandos: List[str] = field(default_factory=list)
    variables: Dict[str, str] = field(default_factory=dict)
    particion: Optional[str] = None
    cuenta: Optional[str] = None
    qos: Optional[str] = None
    nodos: Optional[str] = None
    sbatch_extra: List[str] = field(default_factory=list)
    raiz_proyecto: Optional[str] = None
    email: Optional[str] = None

    # -- construcción -------------------------------------------------------
    @classmethod
    def desde_dict(cls, nombre: str, datos: Dict[str, Any]) -> "PerfilServidor":
        conocidos = {f for f in cls.__dataclass_fields__ if f != "nombre"}
        desconocidos = set(datos) - conocidos
        if desconocidos:
            raise ConfigError(
                f"servidor '{nombre}': claves no reconocidas: "
                + ", ".join(sorted(desconocidos))
            )
        limpio = {k: v for k, v in datos.items() if v is not None or k in ("gres",)}
        obj = cls(nombre=nombre, **limpio)
        obj.validar()
        return obj

    def validar(self) -> None:
        if self.scheduler not in ("local", "slurm"):
            raise ConfigError(
                f"servidor '{self.nombre}': scheduler '{self.scheduler}' no soportado "
                "(usa 'slurm' o 'local')"
            )
        if self.scheduler == "slurm" and not self.particion:
            raise ConfigError(
                f"servidor '{self.nombre}': un perfil de Slurm necesita 'particion'"
            )
        for campo in ("cpus", "memoria_gb", "horas"):
            if int(getattr(self, campo)) <= 0:
                raise ConfigError(f"servidor '{self.nombre}': {campo} debe ser > 0")

    # -- utilidades ---------------------------------------------------------
    @property
    def usa_gpu(self) -> bool:
        return int(self.gpus) > 0

    @property
    def tiempo_slurm(self) -> str:
        return f"{int(self.horas):02d}:00:00"

    def herramientas_requeridas(self) -> List[str]:
        req = ["python"]
        if self.scheduler == "slurm":
            req += ["sbatch", "squeue", "scancel"]
        return req

    def comprobar(self) -> List[Dict[str, Any]]:
        """Comprobaciones de entorno; devuelve una lista de resultados."""
        res: List[Dict[str, Any]] = []
        for h in self.herramientas_requeridas():
            ruta = shutil.which(h)
            res.append(
                {
                    "comprobacion": f"ejecutable '{h}'",
                    "ok": ruta is not None,
                    "detalle": ruta or "no encontrado en PATH",
                }
            )
        if self.conda_env:
            ruta = shutil.which(self.gestor_entorno)
            res.append(
                {
                    "comprobacion": f"gestor de entornos ({self.gestor_entorno})",
                    "ok": ruta is not None,
                    "detalle": ruta or "no encontrado; revisa 'pre_comandos'",
                }
            )
        if self.usa_gpu:
            smi = shutil.which("nvidia-smi")
            detalle = smi or "no encontrado (normal si la GPU sólo existe en el nodo de cómputo)"
            res.append(
                {
                    "comprobacion": "nvidia-smi",
                    "ok": smi is not None,
                    "detalle": detalle,
                    "aviso": self.scheduler == "slurm",
                }
            )
        if self.scheduler == "slurm" and self.particion:
            res.append(self._comprobar_particion())
        return res

    def _comprobar_particion(self) -> Dict[str, Any]:
        if shutil.which("sinfo") is None:
            return {
                "comprobacion": f"partición '{self.particion}'",
                "ok": False,
                "detalle": "sinfo no disponible; no se pudo verificar",
                "aviso": True,
            }
        try:
            salida = subprocess.run(
                ["sinfo", "-h", "-p", str(self.particion), "-o", "%P %c %m %G"],
                capture_output=True,
                text=True,
                timeout=20,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            return {
                "comprobacion": f"partición '{self.particion}'",
                "ok": False,
                "detalle": f"error al consultar sinfo: {exc}",
                "aviso": True,
            }
        texto = (salida.stdout or "").strip()
        if not texto:
            return {
                "comprobacion": f"partición '{self.particion}'",
                "ok": False,
                "detalle": "la partición no existe o no es visible para tu usuario",
            }
        return {
            "comprobacion": f"partición '{self.particion}'",
            "ok": True,
            "detalle": texto.splitlines()[0],
        }


# ---------------------------------------------------------------------------
#  Compuestos
# ---------------------------------------------------------------------------
@dataclass
class Compuesto:
    """Metadatos de un compuesto: lo único que cambia entre simulaciones."""

    clave: str
    nombre: str
    directorio: Path
    #: Donde viven las corridas ya hechas, si no es el propio `directorio`.
    #: En este repositorio las entradas (scripts/, sistema/) se versionan y los
    #: resultados de 50 ns viven aparte, en resultados_50ns/<Sistema>/.
    directorio_resultados: Optional[Path] = None
    pubchem: Optional[int] = None
    titulo: str = ""
    receptor: Optional[str] = None
    ligando: Optional[str] = None
    smiles: Optional[str] = None
    carga_ligando: int = 0
    ligando_resname: str = "LIG"
    receptor_residuos: str = "448-723"
    ph: float = 7.4
    exhaustiveness: int = 100
    caja_vina: List[float] = field(default_factory=lambda: [20.0, 20.0, 20.0])
    bolsillo_residuos: List[int] = field(default_factory=list)
    caja_solvatacion: float = 12.0
    sal_molar: float = 0.15

    @classmethod
    def desde_dict(cls, clave: str, datos: Dict[str, Any], raiz: Path) -> "Compuesto":
        datos = dict(datos)
        directorio = datos.pop("directorio", None)
        if not directorio:
            raise ConfigError(f"compuesto '{clave}': falta 'directorio'")
        ruta = Path(directorio)
        if not ruta.is_absolute():
            ruta = raiz / ruta
        resultados = datos.pop("directorio_resultados", None)
        if resultados:
            resultados = Path(resultados)
            if not resultados.is_absolute():
                resultados = raiz / resultados
            datos["directorio_resultados"] = resultados
        conocidos = {f for f in cls.__dataclass_fields__
                     if f not in ("clave", "directorio")}
        desconocidos = set(datos) - conocidos
        if desconocidos:
            raise ConfigError(
                f"compuesto '{clave}': claves no reconocidas: "
                + ", ".join(sorted(desconocidos))
            )
        datos.setdefault("nombre", clave.capitalize())
        datos.setdefault("titulo", datos["nombre"])
        return cls(clave=clave, directorio=ruta, **datos)

    # -- rutas derivadas ----------------------------------------------------
    @property
    def sistema_dir(self) -> Path:
        return self.directorio / "sistema"

    @property
    def entrada_dir(self) -> Path:
        return self.directorio / "entrada"

    @property
    def runs_dir(self) -> Path:
        return self.base_resultados / "runs"

    @property
    def base_resultados(self) -> Path:
        """Carpeta donde se escriben y se buscan las corridas."""
        return self.directorio_resultados or self.directorio

    @property
    def sistema_preparado(self) -> bool:
        return (self.sistema_dir / "sistema.prmtop").exists()

    def receptor_path(self, raiz: Path) -> Optional[Path]:
        return self._entrada_path(self.receptor, raiz)

    def ligando_path(self, raiz: Path) -> Optional[Path]:
        return self._entrada_path(self.ligando, raiz)

    def _entrada_path(self, nombre: Optional[str], raiz: Path) -> Optional[Path]:
        if not nombre:
            return None
        candidatos = [
            self.entrada_dir / nombre,
            raiz / "simulaciones" / "entrada" / nombre,
            raiz / "Data" / nombre,
            raiz / "data" / nombre,
            raiz.parent / "Data" / nombre,
        ]
        for c in candidatos:
            if c.exists():
                return c
        return candidatos[0]


# ---------------------------------------------------------------------------
#  Proyecto
# ---------------------------------------------------------------------------
@dataclass
class ProyectoConfig:
    """Configuración completa: raíz, servidores, compuestos y tiempos."""

    raiz: Path
    servidores: Dict[str, PerfilServidor]
    compuestos: Dict[str, Compuesto]
    tiempos: Dict[str, Dict[str, Any]]

    # -- acceso -------------------------------------------------------------
    def servidor(self, nombre: str) -> PerfilServidor:
        if nombre in ("auto", "", None):
            nombre = detectar_servidor(self)
        if nombre not in self.servidores:
            raise ConfigError(
                f"servidor '{nombre}' no está en {ARCHIVO_SERVIDORES}. "
                f"Disponibles: {', '.join(sorted(self.servidores))}"
            )
        return self.servidores[nombre]

    def compuesto(self, clave: str) -> Compuesto:
        c = clave.strip().lower()
        if c in self.compuestos:
            return self.compuestos[c]
        # permitir el nombre o el CID de PubChem
        for comp in self.compuestos.values():
            if c in (comp.nombre.lower(), str(comp.pubchem)):
                return comp
        raise ConfigError(
            f"compuesto '{clave}' no está en {ARCHIVO_COMPUESTOS}. "
            f"Disponibles: {', '.join(sorted(self.compuestos))}"
        )

    def resolver_compuestos(self, spec: str) -> List[Compuesto]:
        """``'todos'`` o una lista separada por comas."""
        if spec.strip().lower() in ("todos", "all", "*"):
            return [self.compuestos[k] for k in sorted(self.compuestos)]
        return [self.compuesto(p) for p in spec.split(",") if p.strip()]

    def resolver_tiempos(self, spec: str) -> List[float]:
        """``'10,20,50'``, ``'prueba'`` o ``'todos'`` -> lista de ns."""
        txt = spec.strip().lower()
        if txt in ("todos", "all", "*"):
            return [float(v["ns"]) for _, v in sorted(
                self.tiempos.items(), key=lambda kv: float(kv[1]["ns"]))]
        out: List[float] = []
        for parte in txt.split(","):
            parte = parte.strip()
            if not parte:
                continue
            if parte in self.tiempos:
                out.append(float(self.tiempos[parte]["ns"]))
                continue
            limpio = parte[:-2] if parte.endswith("ns") else parte
            try:
                out.append(float(limpio))
            except ValueError as exc:
                raise ConfigError(
                    f"tiempo '{parte}' no es válido. Usa ns (10, 20, 50, 0.1) "
                    f"o un atajo: {', '.join(sorted(self.tiempos))}"
                ) from exc
        if not out:
            raise ConfigError("no se indicó ningún tiempo de simulación")
        return out

    def parametros_tiempo(self, ns: float) -> Dict[str, Any]:
        """Parámetros derivados del tiempo de producción.

        Si el tiempo coincide con un atajo de ``compuestos.yaml`` se usan sus
        valores; si no, se derivan: segmentos de ~5 ns y un fotograma cada 10 ps.
        """
        for datos in self.tiempos.values():
            if abs(float(datos["ns"]) - float(ns)) < 1e-9:
                base = self._parametros_por_defecto(ns)
                base.update({k: v for k, v in datos.items() if k != "ns"})
                return base
        return self._parametros_por_defecto(ns)

    @staticmethod
    def _parametros_por_defecto(ns: float) -> Dict[str, Any]:
        corto = float(ns) < 1.0
        return {
            "t_produccion_ps": float(ns) * 1000.0,
            "t_calentamiento_ps": 20.0 if corto else 500.0,
            "t_equilibrado_ps": 20.0 if corto else 500.0,
            "segmentos": 4 if corto else max(1, int(round(float(ns) / 5.0))),
            "guardar_cada_ps": 1.0 if corto else 10.0,
            "frames_mmgbsa": 50 if corto else 200,
        }


def cargar_proyecto(raiz: Optional[os.PathLike] = None) -> ProyectoConfig:
    """Lee ``config/servidores.yaml`` y ``config/compuestos.yaml``."""
    base = Path(raiz).expanduser().resolve() if raiz else encontrar_raiz()

    srv_datos = _cargar_yaml(base / ARCHIVO_SERVIDORES)
    srv_defaults = srv_datos.get("defaults") or {}
    servidores: Dict[str, PerfilServidor] = {}
    for nombre, datos in (srv_datos.get("servidores") or {}).items():
        mezcla = dict(srv_defaults)
        mezcla.update(datos or {})
        servidores[nombre] = PerfilServidor.desde_dict(nombre, mezcla)
    if not servidores:
        raise ConfigError(f"{ARCHIVO_SERVIDORES} no define ningún servidor")

    cmp_datos = _cargar_yaml(base / ARCHIVO_COMPUESTOS)
    cmp_defaults = cmp_datos.get("defaults") or {}
    compuestos: Dict[str, Compuesto] = {}
    for clave, datos in (cmp_datos.get("compuestos") or {}).items():
        mezcla = dict(cmp_defaults)
        mezcla.update(datos or {})
        compuestos[clave] = Compuesto.desde_dict(clave, mezcla, base)
    if not compuestos:
        raise ConfigError(f"{ARCHIVO_COMPUESTOS} no define ningún compuesto")

    tiempos = {}
    for clave, datos in (cmp_datos.get("tiempos") or {}).items():
        if not isinstance(datos, dict) or "ns" not in datos:
            raise ConfigError(f"tiempo '{clave}': falta la clave 'ns'")
        tiempos[clave] = datos

    return ProyectoConfig(raiz=base, servidores=servidores,
                          compuestos=compuestos, tiempos=tiempos)


def detectar_servidor(cfg: ProyectoConfig) -> str:
    """Adivina el perfil adecuado para la máquina actual.

    Prioridad: ``$DMLRUN_SERVIDOR``, un perfil cuyo nombre coincida con el
    hostname, el primer perfil de Slurm si ``sbatch`` existe, y ``local``.
    """
    env = os.environ.get("DMLRUN_SERVIDOR")
    if env:
        if env not in cfg.servidores:
            raise ConfigError(
                f"$DMLRUN_SERVIDOR={env} no está en {ARCHIVO_SERVIDORES}"
            )
        return env

    host = socket.gethostname().split(".")[0].lower()
    for nombre in cfg.servidores:
        if nombre.lower() == host:
            return nombre

    if shutil.which("sbatch"):
        for nombre, perfil in cfg.servidores.items():
            if perfil.scheduler == "slurm":
                return nombre
    # Sin cola a la vista, el perfil por defecto debe ser uno que ejecute aqui
    # mismo: devolver el primero del archivo podria ser un perfil de Slurm y
    # haria fallar el envio en una maquina que no tiene sbatch.
    for preferido in ("local", "actual", "portatil"):
        if preferido in cfg.servidores:
            return preferido
    for nombre, perfil in cfg.servidores.items():
        if perfil.scheduler == "local":
            return nombre
    return next(iter(cfg.servidores))
