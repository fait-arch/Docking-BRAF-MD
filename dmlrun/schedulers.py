# -*- coding: utf-8 -*-
"""
Backends de ejecución: Slurm y ejecución directa.

Ambos reciben el mismo :class:`JobSpec` (comandos + recursos + entorno) y
generan un script de shell equivalente; lo único que cambia es la cabecera y
la forma de lanzarlo:

* :class:`SlurmScheduler` escribe las directivas ``#SBATCH`` y llama a ``sbatch``.
* :class:`LocalScheduler` lanza el script con ``setsid`` y guarda el PID.

En los dos casos el script envuelve los comandos con un ``trap`` que escribe el
código de salida en ``estado/codigo_salida``, de modo que ``dmlrun status``
funciona igual con o sin gestor de colas.
"""
from __future__ import annotations

import os
import shlex
import subprocess
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional

ESTADOS_ACTIVOS = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED",
                   "REQUEUED", "RESIZING"}


class SchedulerError(RuntimeError):
    """Fallo al enviar o consultar un trabajo."""


@dataclass
class JobSpec:
    """Descripción de un trabajo, independiente del gestor de colas."""

    nombre: str
    comandos: List[str]
    workdir: Path
    log: Path
    script: Path
    estado_dir: Path
    cpus: int = 4
    memoria_gb: int = 16
    horas: int = 8
    gpus: int = 0
    gres: Optional[str] = None
    particion: Optional[str] = None
    cuenta: Optional[str] = None
    qos: Optional[str] = None
    nodos: Optional[str] = None
    email: Optional[str] = None
    modulos: List[str] = field(default_factory=list)
    conda_env: Optional[str] = None
    gestor_entorno: str = "micromamba"
    pre_comandos: List[str] = field(default_factory=list)
    post_comandos: List[str] = field(default_factory=list)
    variables: Dict[str, str] = field(default_factory=dict)
    sbatch_extra: List[str] = field(default_factory=list)

    def cuerpo(self) -> List[str]:
        """Las líneas comunes a los dos backends (entorno + comandos)."""
        out: List[str] = ["set -euo pipefail", ""]
        out += [
            f'ESTADO={shlex.quote(str(self.estado_dir))}',
            'mkdir -p "$ESTADO"',
            'date -Is > "$ESTADO/inicio"',
            'trap \'c=$?; echo "$c" > "$ESTADO/codigo_salida"; '
            'date -Is > "$ESTADO/fin"; exit $c\' EXIT',
            "",
            'echo ""',
            'echo "==============================================================="',
            "echo \"== servidor  : $(hostname)\"",
            "echo \"== fecha     : $(date -Is)\"",
            "",
        ]
        for m in self.modulos:
            out.append(f"module load {m}")
        if self.modulos:
            out.append("")
        for c in self.pre_comandos:
            out.append(c)
        if self.conda_env:
            out.append(f"{self.gestor_entorno} activate {shlex.quote(self.conda_env)}")
        if self.pre_comandos or self.conda_env:
            out.append("")
        for k, v in self.variables.items():
            out.append(f'export {k}="{v}"')
        if self.variables:
            out.append("")
        out.append(f"cd {shlex.quote(str(self.workdir))}")
        out.append('echo "== python    : $(command -v python)"')
        out.append("")
        out += list(self.comandos)
        if self.post_comandos:
            out.append("")
            out += list(self.post_comandos)
        out.append("")
        out.append('echo "== terminado : $(date -Is)"')
        return out


@dataclass
class SubmitResult:
    """Resultado de enviar un trabajo."""

    job_id: str
    scheduler: str
    script: Path
    log: Path
    detalle: str = ""


class Scheduler:
    """Interfaz común."""

    nombre = "base"

    def render(self, job: JobSpec) -> str:  # pragma: no cover
        raise NotImplementedError

    def escribir_script(self, job: JobSpec) -> Path:
        job.script.parent.mkdir(parents=True, exist_ok=True)
        texto = self.render(job)
        job.script.write_text(texto, encoding="utf-8")
        job.script.chmod(0o755)
        return job.script

    def submit(self, job: JobSpec) -> SubmitResult:  # pragma: no cover
        raise NotImplementedError

    def status(self, job_id: str, estado_dir: Optional[Path] = None) -> str:  # pragma: no cover
        raise NotImplementedError

    def cancel(self, job_id: str) -> str:  # pragma: no cover
        raise NotImplementedError

    # -- común --------------------------------------------------------------
    @staticmethod
    def _estado_por_archivos(estado_dir: Optional[Path]) -> Optional[str]:
        """Estado deducido de los archivos que escribe el propio trabajo.

        Devuelve ``None`` mientras el trabajo siga en marcha. El ``codigo_salida``
        de una ejecución anterior no se confunde con el de la actual: se compara
        con la marca ``inicio``, que el script reescribe al arrancar (borrarlo no
        es una opción fiable, porque hay carpetas donde no se puede borrar).
        """
        if not estado_dir:
            return None
        base = Path(estado_dir)
        salida, inicio = base / "codigo_salida", base / "inicio"
        if not salida.exists():
            return None
        try:
            if inicio.exists() and inicio.stat().st_mtime > salida.stat().st_mtime:
                return None          # el trabajo actual todavía no ha terminado
        except OSError:
            pass
        try:
            codigo = int(salida.read_text().strip() or "1")
        except (OSError, ValueError):
            return "DESCONOCIDO"
        return "COMPLETADO" if codigo == 0 else f"FALLIDO (código {codigo})"


class SlurmScheduler(Scheduler):
    """Envía trabajos con ``sbatch`` y los consulta con ``squeue``/``sacct``."""

    nombre = "slurm"

    def render(self, job: JobSpec) -> str:
        h: List[str] = ["#!/bin/bash"]
        h.append(f"#SBATCH --job-name={job.nombre}")
        if job.particion:
            h.append(f"#SBATCH --partition={job.particion}")
        if job.cuenta:
            h.append(f"#SBATCH --account={job.cuenta}")
        if job.qos:
            h.append(f"#SBATCH --qos={job.qos}")
        if job.nodos:
            h.append(f"#SBATCH --nodelist={job.nodos}")
        h.append("#SBATCH --nodes=1")
        h.append("#SBATCH --ntasks=1")
        h.append(f"#SBATCH --cpus-per-task={int(job.cpus)}")
        h.append(f"#SBATCH --mem={int(job.memoria_gb)}G")
        h.append(f"#SBATCH --time={int(job.horas):02d}:00:00")
        if int(job.gpus) > 0:
            # Si el perfil indica un GRES explícito se usa tal cual; es lo que
            # distingue 'gpu:1' de un perfil MIG como 'gpu:2g.10gb'.
            h.append(f"#SBATCH --gres={job.gres or f'gpu:{int(job.gpus)}'}")
        h.append(f"#SBATCH --output={job.log}")
        h.append(f"#SBATCH --error={job.log}")
        if job.email:
            h.append(f"#SBATCH --mail-user={job.email}")
            h.append("#SBATCH --mail-type=END,FAIL")
        for extra in job.sbatch_extra:
            linea = str(extra).strip()
            if not linea:
                continue
            h.append(linea if linea.startswith("#SBATCH") else f"#SBATCH {linea}")
        h.append("")
        return "\n".join(h + job.cuerpo()) + "\n"

    def submit(self, job: JobSpec) -> SubmitResult:
        script = self.escribir_script(job)
        job.log.parent.mkdir(parents=True, exist_ok=True)
        try:
            proc = subprocess.run(
                ["sbatch", "--parsable", str(script)],
                capture_output=True, text=True, cwd=str(job.workdir), timeout=120,
            )
        except FileNotFoundError as exc:
            raise SchedulerError(
                "No encuentro 'sbatch'. ¿Estás en el nodo de login de la HPC? "
                "Para un servidor sin cola usa un perfil con scheduler: local."
            ) from exc
        if proc.returncode != 0:
            raise SchedulerError(self._explicar_error(proc.stderr or proc.stdout, job))
        job_id = (proc.stdout or "").strip().split(";")[0]
        if not job_id:
            raise SchedulerError(f"sbatch no devolvió un id de trabajo: {proc.stdout!r}")
        return SubmitResult(job_id=job_id, scheduler=self.nombre,
                            script=script, log=job.log)

    @staticmethod
    def _explicar_error(stderr: str, job: JobSpec) -> str:
        txt = (stderr or "").strip()
        msg = [f"sbatch rechazó el trabajo:\n  {txt}"]
        if "Requested node configuration is not available" in txt:
            msg.append(
                "\nEsto significa que ningún nodo de la partición existe con esa "
                "combinación. Revisa, en este orden:\n"
                f"  1. GPU: pediste gres='{job.gres or f'gpu:{job.gpus}'}'. Comprueba el nombre real con\n"
                f"       scontrol show node <nodo> | grep -i gres\n"
                f"  2. Recursos: {job.cpus} CPU y {job.memoria_gb} GB. Comprueba el máximo con\n"
                f"       sinfo -p {job.particion} -o '%P %c %m %G'\n"
                "  3. Tiempo: el límite de la partición puede ser menor que --horas.\n"
                "Ajusta el perfil en config/servidores.yaml y vuelve a enviar."
            )
        elif "Invalid partition" in txt or "invalid partition" in txt:
            msg.append(
                f"\nLa partición '{job.particion}' no existe o no está disponible para tu "
                "usuario. Lista las válidas con: sinfo -s"
            )
        elif "Invalid account" in txt:
            msg.append("\nRevisa 'cuenta' en el perfil del servidor (sacctmgr show assoc user=$USER).")
        return "\n".join(msg)

    def status(self, job_id: str, estado_dir: Optional[Path] = None) -> str:
        try:
            proc = subprocess.run(["squeue", "-h", "-j", str(job_id), "-o", "%T"],
                                  capture_output=True, text=True, timeout=30)
            estado = (proc.stdout or "").strip().splitlines()
            if estado and estado[0].strip():
                return estado[0].strip()
        except (OSError, subprocess.SubprocessError):
            pass
        try:
            proc = subprocess.run(
                ["sacct", "-n", "-X", "-j", str(job_id), "-o", "State"],
                capture_output=True, text=True, timeout=30)
            estado = (proc.stdout or "").strip().splitlines()
            if estado and estado[0].strip():
                return estado[0].strip().split()[0]
        except (OSError, subprocess.SubprocessError):
            pass
        return self._estado_por_archivos(estado_dir) or "DESCONOCIDO"

    def cancel(self, job_id: str) -> str:
        proc = subprocess.run(["scancel", str(job_id)], capture_output=True, text=True)
        if proc.returncode != 0:
            raise SchedulerError(f"scancel falló: {proc.stderr.strip()}")
        return f"trabajo {job_id} cancelado"


class LocalScheduler(Scheduler):
    """Ejecuta en el propio servidor, en segundo plano y desacoplado del shell."""

    nombre = "local"

    def render(self, job: JobSpec) -> str:
        h = ["#!/bin/bash",
             f"# trabajo: {job.nombre}",
             f"# generado: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
             ""]
        return "\n".join(h + job.cuerpo()) + "\n"

    def submit(self, job: JobSpec) -> SubmitResult:
        script = self.escribir_script(job)
        job.log.parent.mkdir(parents=True, exist_ok=True)
        pid_file = job.estado_dir / "pid"
        job.estado_dir.mkdir(parents=True, exist_ok=True)
        # setsid desacopla el proceso de la terminal: sigue vivo al cerrar la sesión.
        with open(job.log, "ab") as log:
            proc = subprocess.Popen(
                ["setsid", "bash", str(script)],
                stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                cwd=str(job.workdir), start_new_session=True,
            )
        pid_file.write_text(str(proc.pid), encoding="utf-8")
        return SubmitResult(job_id=str(proc.pid), scheduler=self.nombre,
                            script=script, log=job.log,
                            detalle="ejecución directa en segundo plano")

    def status(self, job_id: str, estado_dir: Optional[Path] = None) -> str:
        por_archivo = self._estado_por_archivos(estado_dir)
        if por_archivo:
            return por_archivo
        try:
            pid = int(job_id)
        except (TypeError, ValueError):
            return "DESCONOCIDO"
        if Path(f"/proc/{pid}").exists():
            return "RUNNING"
        try:
            os.kill(pid, 0)
            return "RUNNING"
        except ProcessLookupError:
            return "INTERRUMPIDO"
        except PermissionError:
            return "RUNNING"
        except OSError:
            return "DESCONOCIDO"

    def cancel(self, job_id: str) -> str:
        import signal
        try:
            pid = int(job_id)
        except (TypeError, ValueError) as exc:
            raise SchedulerError(f"PID inválido: {job_id}") from exc
        try:
            os.killpg(os.getpgid(pid), signal.SIGTERM)
        except ProcessLookupError as exc:
            raise SchedulerError(f"el proceso {pid} ya no existe") from exc
        except OSError:
            os.kill(pid, signal.SIGTERM)
        return f"proceso {pid} terminado"


def obtener_scheduler(nombre: str) -> Scheduler:
    """``'slurm'`` o ``'local'`` -> instancia del backend."""
    if nombre == "slurm":
        return SlurmScheduler()
    if nombre == "local":
        return LocalScheduler()
    raise SchedulerError(f"scheduler '{nombre}' no soportado")
