# -*- coding: utf-8 -*-
"""
CLI de orquestación: ``dmlrun``.

    dmlrun servidores                      lista los perfiles de entorno
    dmlrun compuestos                      lista los compuestos registrados
    dmlrun config --servidor hpc           muestra y comprueba un perfil
    dmlrun run     ...                     envía ejecuciones
    dmlrun status  ...                     consulta el estado
    dmlrun logs    --id ID                 muestra el log
    dmlrun cancel  --id ID                 cancela
    dmlrun comparar --compuesto pallidine  compara tiempos de simulación
    dmlrun reporte  --compuesto todos      informe HTML con figuras y datos

Lo único obligatorio para cambiar de entorno es ``--servidor``: el resto
(scheduler, recursos, módulos, plataforma de OpenMM) sale de
``config/servidores.yaml``.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .config import ConfigError, ProyectoConfig, cargar_proyecto, detectar_servidor
from .layout import NS_PRUEBA, RunLayout, etiqueta_tiempo
from .pipelines import PIPELINES, Pipeline, construir_pipeline
from .registry import (Registro, RunRecord, construir_manifiesto,
                       escribir_manifiesto, huella_parametros)
from .schedulers import JobSpec, SchedulerError, obtener_scheduler

PIPELINE_ALIAS = {
    "resultados": ["resultados", "results", "pipelineresults", "md"],
    "ml": ["ml", "machinelearning", "machinelearningpipeline", "rescoring"],
}


# ---------------------------------------------------------------------------
#  Presentación
# ---------------------------------------------------------------------------
def _tabla(filas: List[Sequence[str]], cabecera: Sequence[str]) -> str:
    if not filas:
        return "  (vacío)"
    cols = len(cabecera)
    ancho = [len(str(cabecera[i])) for i in range(cols)]
    for f in filas:
        for i in range(cols):
            ancho[i] = max(ancho[i], len(str(f[i]) if i < len(f) else ""))
    sep = "  "
    out = [sep.join(str(cabecera[i]).ljust(ancho[i]) for i in range(cols)).rstrip()]
    out.append(sep.join("-" * ancho[i] for i in range(cols)))
    for f in filas:
        out.append(sep.join(
            str(f[i] if i < len(f) else "").ljust(ancho[i]) for i in range(cols)
        ).rstrip())
    return "\n".join("  " + l for l in out)


def _marca(ok: bool, aviso: bool = False) -> str:
    if ok:
        return "ok"
    return "aviso" if aviso else "FALTA"


def _resolver_pipeline(valor: str) -> List[str]:
    v = valor.strip().lower()
    if v in ("ambos", "both", "todos", "all"):
        return ["resultados", "ml"]
    for clave, alias in PIPELINE_ALIAS.items():
        if v in alias:
            return [clave]
    raise ConfigError(
        f"pipeline '{valor}' desconocido. Usa 'resultados', 'ml' o 'ambos'"
    )


# ---------------------------------------------------------------------------
#  Subcomandos informativos
# ---------------------------------------------------------------------------
def cmd_servidores(cfg: ProyectoConfig, args) -> int:
    filas = []
    actual = detectar_servidor(cfg)
    for nombre in sorted(cfg.servidores):
        p = cfg.servidores[nombre]
        filas.append([
            nombre + (" *" if nombre == actual else ""),
            p.scheduler,
            p.particion or "-",
            str(p.cpus),
            f"{p.memoria_gb} GB",
            f"{p.gpus}" + (f" ({p.gres})" if p.gres and p.gpus else ""),
            p.plataforma_openmm,
            p.descripcion,
        ])
    print("Perfiles de servidor (config/servidores.yaml):\n")
    print(_tabla(filas, ["servidor", "cola", "partición", "cpu", "memoria",
                         "gpu", "openmm", "descripción"]))
    print(f"\n  * perfil que se usaría con --servidor auto en esta máquina")
    return 0


def cmd_compuestos(cfg: ProyectoConfig, args) -> int:
    filas = []
    for clave in sorted(cfg.compuestos):
        c = cfg.compuestos[clave]
        runs = []
        if c.runs_dir.exists():
            runs = sorted(d.name for d in c.runs_dir.iterdir() if d.is_dir())
        # tiempos en carpetas heredadas
        for d in sorted(c.base_resultados.glob("md_*")):
            etq = d.name[3:]
            etq = "0.1ns" if etq == "prueba" else etq
            if etq not in runs:
                runs.append(etq + " (heredado)")
        filas.append([
            clave,
            c.nombre,
            str(c.pubchem or "-"),
            "sí" if c.sistema_preparado else "no",
            ", ".join(runs) or "-",
        ])
    print("Compuestos registrados (config/compuestos.yaml):\n")
    print(_tabla(filas, ["clave", "nombre", "pubchem", "sistema", "ejecuciones"]))
    print("\n  Tiempos predefinidos: " +
          ", ".join(f"{k} ({v['ns']:g} ns)" for k, v in sorted(
              cfg.tiempos.items(), key=lambda kv: float(kv[1]["ns"]))))
    return 0


def cmd_config(cfg: ProyectoConfig, args) -> int:
    perfil = cfg.servidor(args.servidor)
    if args.json:
        from dataclasses import asdict
        print(json.dumps({"raiz": str(cfg.raiz), "servidor": asdict(perfil)},
                         indent=2, ensure_ascii=False))
        return 0

    print(f"Raíz del proyecto : {cfg.raiz}")
    print(f"Servidor          : {perfil.nombre}  —  {perfil.descripcion}\n")
    filas = [
        ["gestor de colas", perfil.scheduler],
        ["partición", perfil.particion or "-"],
        ["cpus", str(perfil.cpus)],
        ["memoria", f"{perfil.memoria_gb} GB"],
        ["tiempo máximo", f"{perfil.horas} h"],
        ["gpus", f"{perfil.gpus}" + (f"  (gres: {perfil.gres})" if perfil.gres else "")],
        ["plataforma OpenMM", perfil.plataforma_openmm],
        ["entorno conda", perfil.conda_env or "-"],
        ["módulos", ", ".join(perfil.modulos) or "-"],
        ["variables", ", ".join(f"{k}={v}" for k, v in perfil.variables.items()) or "-"],
    ]
    print(_tabla(filas, ["parámetro", "valor"]))

    print("\nComprobaciones en esta máquina:\n")
    res = perfil.comprobar()
    print(_tabla([[r["comprobacion"], _marca(r["ok"], r.get("aviso", False)),
                   r["detalle"]] for r in res],
                 ["comprobación", "estado", "detalle"]))
    criticos = [r for r in res if not r["ok"] and not r.get("aviso")]
    if criticos:
        print("\n  Faltan herramientas críticas para este perfil. Si estás en otra "
              "máquina, usa\n  el perfil que corresponda o ejecuta con --dry-run "
              "para generar sólo el script.")
        return 1
    return 0


# ---------------------------------------------------------------------------
#  run
# ---------------------------------------------------------------------------
#: Carpeta de la raíz del repo donde se crean las corridas nuevas.
DIR_ENTREGAS = "resultados"

#: Sello de la sesión: se calcula UNA vez por invocación de la CLI, para que
#: todos los compuestos lanzados juntos caigan en la misma carpeta de fecha.
#: Si se calculara por compuesto, una tanda que cruce un minuto se partiría.
SELLO_SESION = datetime.now().strftime("%Y%m%d-%H%M")

#: Una carpeta de sesión: AAAAMMDD-HHMM.
RE_SESION = re.compile(r"^\d{8}-\d{4}$")


def _nombre_entrega(compuesto, etiqueta: str) -> str:
    """``'10ns'`` + ``'Notholaenic acid'`` -> ``'10nsNotholaenicAcid'``.

    Sin espacios a propósito: las directivas ``#SBATCH`` se parten por espacios
    y no admiten comillas, así que un ``--output=.../10nsNotholaenic acid/...``
    rompe el envío a Slurm.
    """
    return etiqueta + "".join(p[:1].upper() + p[1:] for p in compuesto.nombre.split())


def _carpeta_entrega(cfg, compuesto, etiqueta: str) -> Path:
    """``resultados/<sesión>/<tiempo>ns<Nombre>``."""
    return (cfg.raiz / DIR_ENTREGAS / SELLO_SESION /
            _nombre_entrega(compuesto, etiqueta))


def _carpetas_de_compuesto(cfg) -> List[Path]:
    """Todas las carpetas de corrida bajo resultados/, sean de la sesión que
    sean, incluidas las heredadas que cuelgan directamente de la raíz."""
    base = cfg.raiz / DIR_ENTREGAS
    if not base.is_dir():
        return []
    salida = []
    for d in base.iterdir():
        if not d.is_dir():
            continue
        if RE_SESION.match(d.name):
            salida += [x for x in d.iterdir() if x.is_dir()]
        else:
            salida.append(d)
    return salida


def _corridas_previas(cfg, compuesto, etiqueta: str) -> List[Path]:
    """Corridas ya existentes de ese compuesto y tiempo, de la más nueva a la
    más vieja.

    Reconoce las tres disposiciones que ha tenido el proyecto: la actual
    (``<sesión>/<tiempo>ns<Nombre>``), la intermedia con la fecha pegada al
    nombre, y la plana original. Así ninguna corrida vieja queda huérfana.
    """
    nombre = _nombre_entrega(compuesto, etiqueta)
    hallados = []
    for d in _carpetas_de_compuesto(cfg):
        if d.name == nombre:
            sello = d.parent.name if RE_SESION.match(d.parent.name) else "00000000-0000"
            hallados.append((sello, d))
            continue
        m = re.fullmatch(rf"{re.escape(nombre)}_(\d{{8}}-\d{{4}})", d.name)
        if m:
            hallados.append((m.group(1), d))
    return [c for _, c in sorted(hallados, key=lambda x: x[0], reverse=True)]


def _topologia_rehecha(compuesto, carpeta: Path) -> bool:
    """¿Se reconstruyó el sistema DESPUÉS de empezar esa corrida?

    Si es así, reanudar sería mezclar una trayectoria vieja con una topología
    nueva: el pipeline daría por buenos los .xml existentes y el resultado
    sería un híbrido silencioso. En ese caso se abre una carpeta nueva.
    """
    try:
        prmtop = (compuesto.sistema_dir / "sistema.prmtop").stat().st_mtime
        marca = (carpeta / "md" / "01_min.xml").stat().st_mtime
    except OSError:
        return False
    return prmtop > marca


def _layout(cfg, compuesto, ns: float, modo: str = "auto",
            nueva: bool = False) -> RunLayout:
    """Layout de una corrida.

    Las corridas nuevas viven en una carpeta plana propia dentro de
    ``resultados/`` -- ``resultados/10nsMacaridine/{md,analisis,resultados,
    ml,estado}`` -- autocontenida y lista para copiar tal cual.

    Lo ya calculado no se mueve: si un tiempo tiene datos en el árbol canónico
    del compuesto (``resultados_50ns/<Sistema>/md_50ns`` y compañía, de las
    corridas heredadas), se sigue leyendo de ahí. ``--layout runs`` o
    ``legacy`` fuerzan el comportamiento anterior.
    """
    etiqueta = etiqueta_tiempo(ns)
    canonico = compuesto.base_resultados
    if modo in ("runs", "legacy"):
        return RunLayout.resolver(canonico, ns, modo=modo,
                                  entradas=compuesto.directorio)
    # Se reanuda la corrida más reciente salvo que se pida una nueva, o que el
    # sistema se haya reconstruido desde entonces. Nada se sobrescribe jamás:
    # cada ejecución nueva estrena su propia carpeta fechada.
    previas = _corridas_previas(cfg, compuesto, etiqueta)
    reutilizable = next((c for c in previas
                         if not _topologia_rehecha(compuesto, c)), None)
    plana = (_carpeta_entrega(cfg, compuesto, etiqueta)
             if (nueva or reutilizable is None) else reutilizable)
    ya_en_canonico = ((canonico / f"md_{etiqueta}").exists()
                      or (canonico / "runs" / etiqueta / "md").exists()
                      or (abs(float(ns) - NS_PRUEBA) < 1e-9
                          and (canonico / "md_prueba").exists()))
    if not previas and ya_en_canonico:
        return RunLayout.resolver(canonico, ns, modo="auto",
                                  entradas=compuesto.directorio)
    return RunLayout(
        base=plana, ns=float(ns), etiqueta=etiqueta,
        md=plana / "md", analisis=plana / "analisis",
        resultados=plana / "resultados", ml=plana / "ml",
        estado=plana / "estado", legacy=False,
        entradas=compuesto.directorio, plano=True,
    )


def _opciones_pipeline(args) -> Dict:
    return {
        "t_calentamiento_ps": args.t_calentamiento,
        "t_equilibrado_ps": args.t_equilibrado,
        "segmentos": args.segmentos,
        "guardar_cada_ps": args.guardar_cada,
        "frames_mmgbsa": args.frames_mmgbsa,
        "max_segmentos": args.max_segmentos,
        "semilla": args.semilla,
        "fuente": args.fuente,
        "motor": args.motor,
        "modelo": args.modelo,
        "frames_ml": args.frames_ml,
    }


def _job_de_pipeline(pipe: Pipeline, perfil, layout: RunLayout,
                     etapas: Optional[List[str]], forzar: bool,
                     horas: Optional[int], cpus: Optional[int],
                     memoria: Optional[int]) -> JobSpec:
    comandos = pipe.comandos(etapas, forzar=forzar)
    return JobSpec(
        nombre=pipe.etiqueta_trabajo(),
        comandos=comandos,
        workdir=pipe.compuesto.directorio,
        log=layout.log,
        script=layout.estado / f"trabajo-{pipe.clave}.sh",
        estado_dir=layout.estado,
        cpus=int(cpus or perfil.cpus),
        memoria_gb=int(memoria or perfil.memoria_gb),
        horas=int(horas or perfil.horas),
        gpus=int(perfil.gpus),
        gres=perfil.gres,
        particion=perfil.particion,
        cuenta=perfil.cuenta,
        qos=perfil.qos,
        nodos=perfil.nodos,
        email=perfil.email,
        modulos=list(perfil.modulos),
        conda_env=perfil.conda_env,
        gestor_entorno=perfil.gestor_entorno,
        pre_comandos=list(perfil.pre_comandos),
        post_comandos=list(perfil.post_comandos),
        variables=dict(perfil.variables),
        sbatch_extra=list(perfil.sbatch_extra),
    )


def cmd_run(cfg: ProyectoConfig, args) -> int:
    perfil = cfg.servidor(args.servidor)
    pipelines = _resolver_pipeline(args.pipeline)
    compuestos = cfg.resolver_compuestos(args.compuesto)
    tiempos = cfg.resolver_tiempos(args.tiempo)
    etapas = [e.strip() for e in args.etapas.split(",")] if args.etapas else None
    registro = Registro(cfg.raiz)
    scheduler = obtener_scheduler(perfil.scheduler)
    opciones = _opciones_pipeline(args)

    print(f"Servidor : {perfil.nombre} ({perfil.scheduler}"
          + (f", partición {perfil.particion}" if perfil.particion else "") + ")")
    print(f"Pipelines: {', '.join(pipelines)}")
    print(f"Compuestos: {', '.join(c.nombre for c in compuestos)}")
    print(f"Tiempos  : {', '.join(etiqueta_tiempo(t) for t in tiempos)}")
    if args.dry_run:
        print("Modo      : --dry-run (no se envía nada)")
    print()

    enviados: List[RunRecord] = []
    problemas = 0

    for compuesto in compuestos:
        for ns in tiempos:
            layout = _layout(cfg, compuesto, ns, modo=args.layout,
                                 nueva=args.nueva_corrida)
            for clave in pipelines:
                pipe = construir_pipeline(clave, cfg, compuesto, layout, perfil,
                                          opciones)
                cabecera = (f"[{compuesto.clave} · {clave} · {layout.etiqueta}]"
                            + (" (carpetas heredadas)" if layout.legacy else ""))
                print(cabecera)

                faltan = getattr(pipe, "requisitos", lambda: [])()
                if faltan:
                    for f in faltan:
                        print(f"  requisito: {f}")
                    if not args.ignorar_requisitos:
                        print("  omitido (usa --ignorar-requisitos para enviarlo "
                              "igualmente)\n")
                        problemas += 1
                        continue

                try:
                    plan = pipe.plan(etapas, forzar=args.forzar)
                except ValueError as exc:
                    print(f"  error: {exc}\n")
                    return 2
                for e in plan:
                    estado = f"omitida ({e.omitida_por})" if e.omitida else "pendiente"
                    print(f"    {e.nombre:<9} {estado:<46} {e.descripcion}")
                if all(e.omitida for e in plan):
                    print("  nada que hacer; usa --forzar para repetir\n")
                    continue

                if not args.dry_run:
                    # El directorio del compuesto es el workdir del trabajo: para un
                    # compuesto recien registrado todavia no existe, y tanto el `cd`
                    # del script como el cwd= de Popen fallarian antes de empezar.
                    compuesto.directorio.mkdir(parents=True, exist_ok=True)
                    layout.crear()
                job = _job_de_pipeline(pipe, perfil, layout, etapas, args.forzar,
                                       args.horas, args.cpus, args.memoria)

                if args.dry_run:
                    script = scheduler.render(job)
                    if args.mostrar_script:
                        print("\n" + "-" * 72)
                        print(script)
                        print("-" * 72)
                    else:
                        destino = layout.estado / f"trabajo-{clave}.sh"
                        print(f"  script que se habría enviado: {destino}")
                        print("  (añade --mostrar-script para verlo aquí)")
                    print()
                    continue

                manifiesto = construir_manifiesto(
                    raiz=cfg.raiz, compuesto=compuesto, layout=layout,
                    servidor=perfil, pipeline=clave,
                    etapas=[e.nombre for e in plan if not e.omitida],
                    parametros=(pipe.params if clave == "resultados"
                                else {"fuente": pipe.fuente, "motor": pipe.motor,
                                      "frames_ml": opciones["frames_ml"]}),
                )

                if args.solo_script:
                    ruta = scheduler.escribir_script(job)
                    escribir_manifiesto(layout, manifiesto)
                    print(f"  script escrito (no enviado): {ruta}")
                    cmd = (f"sbatch {ruta}" if perfil.scheduler == "slurm"
                           else f"setsid bash {ruta} &")
                    print(f"  envíalo cuando quieras con:  {cmd}\n")
                    continue

                try:
                    res = scheduler.submit(job)
                except SchedulerError as exc:
                    print(f"  ERROR al enviar:\n{_sangrar(str(exc))}\n")
                    problemas += 1
                    continue

                run_id = (f"{compuesto.clave}-{clave}-{layout.etiqueta}-"
                          f"{datetime.now().strftime('%Y%m%d-%H%M%S')}")
                record = RunRecord(
                    id=run_id,
                    compuesto=compuesto.clave,
                    pipeline=clave,
                    tiempo_ns=float(ns),
                    etiqueta=layout.etiqueta,
                    servidor=perfil.nombre,
                    scheduler=perfil.scheduler,
                    job_id=res.job_id,
                    estado="ENVIADO",
                    creado=datetime.now().astimezone().isoformat(timespec="seconds"),
                    directorio=str(layout.raiz_run),
                    log=str(res.log),
                    script=str(res.script),
                    estado_dir=str(layout.estado),
                    etapas=[e.nombre for e in plan if not e.omitida],
                    huella=manifiesto["huella"],
                    git_commit=manifiesto["codigo"]["git_commit"],
                )
                manifiesto["run_id"] = run_id
                manifiesto["job_id"] = res.job_id
                escribir_manifiesto(layout, manifiesto)
                registro.añadir(record)
                enviados.append(record)
                etiqueta_id = "job" if perfil.scheduler == "slurm" else "pid"
                print(f"  enviado: {etiqueta_id} {res.job_id}  |  id {run_id}")
                print(f"  log    : {res.log}\n")

    if enviados:
        print(f"{len(enviados)} ejecución(es) enviada(s). Consulta el estado con:")
        print("    dmlrun status")
    if problemas:
        print(f"\n{problemas} ejecución(es) no se enviaron (ver arriba).")
        return 1
    return 0


def _sangrar(texto: str, prefijo: str = "    ") -> str:
    return "\n".join(prefijo + l for l in texto.splitlines())


# ---------------------------------------------------------------------------
#  status / logs / cancel
# ---------------------------------------------------------------------------
def cmd_status(cfg: ProyectoConfig, args) -> int:
    registro = Registro(cfg.raiz)
    registros = registro.cargar()
    if args.id:
        registros = [registro.buscar(args.id)]
    if args.compuesto:
        claves = {c.clave for c in cfg.resolver_compuestos(args.compuesto)}
        registros = [r for r in registros if r.compuesto in claves]
    if args.pipeline:
        claves = set(_resolver_pipeline(args.pipeline))
        registros = [r for r in registros if r.pipeline in claves]
    if not registros:
        print("No hay ejecuciones registradas todavía.")
        print("Lanza una con:  dmlrun run --pipeline resultados --compuesto "
              "pallidine --tiempo 10 --servidor hpc")
        return 0

    if not args.todas:
        registros = registros[-args.n:]

    filas = []
    for r in registros:
        sched = obtener_scheduler(r.scheduler)
        estado = sched.status(r.job_id, Path(r.estado_dir))
        if estado != r.estado:
            registro.actualizar_estado(r.id, estado)
        progreso = _progreso(cfg, r)
        filas.append([r.id, r.compuesto, r.pipeline, r.etiqueta, r.servidor,
                      r.job_id, estado, progreso])
    print(_tabla(filas, ["id", "compuesto", "pipeline", "tiempo", "servidor",
                         "job/pid", "estado", "progreso"]))

    if args.id and len(registros) == 1:
        r = registros[0]
        print(f"\n  directorio : {r.directorio}")
        print(f"  log        : {r.log}")
        print(f"  script     : {r.script}")
        print(f"  etapas     : {', '.join(r.etapas) or '-'}")
        print(f"  huella     : {r.huella}  (parámetros científicos)")
        print(f"  commit     : {r.git_commit or '-'}")
    return 0


def _progreso(cfg: ProyectoConfig, r: RunRecord) -> str:
    """Progreso real leído del disco, no del gestor de colas."""
    try:
        compuesto = cfg.compuesto(r.compuesto)
    except ConfigError:
        return "-"
    layout = _layout(cfg, compuesto, r.tiempo_ns)
    if r.pipeline == "ml":
        hitos = [("features", layout.ml / "features.csv"),
                 ("rescored", layout.ml / "rescored.csv")]
        hechos = [n for n, p in hitos if p.exists()]
        return ", ".join(hechos) or "sin salida aún"
    total = int(cfg.parametros_tiempo(r.tiempo_ns)["segmentos"])
    hechos = layout.segmentos_hechos()
    etapas = layout.etapas_completadas()
    if hechos:
        return f"producción {hechos}/{total}"
    for nombre in ("eq", "cal", "min"):
        if etapas.get(nombre):
            return f"{nombre} hecho"
    return "iniciando"


def cmd_logs(cfg: ProyectoConfig, args) -> int:
    registro = Registro(cfg.raiz)
    r = registro.buscar(args.id)
    log = Path(r.log)
    if not log.exists():
        print(f"El log todavía no existe: {log}")
        return 1
    if args.seguir:
        try:
            subprocess.run(["tail", "-n", str(args.n), "-f", str(log)])
        except KeyboardInterrupt:
            pass
        return 0
    lineas = log.read_text(errors="replace").splitlines()
    print(f"# {log}  (últimas {min(args.n, len(lineas))} de {len(lineas)} líneas)\n")
    print("\n".join(lineas[-args.n:]))
    return 0


def cmd_cancel(cfg: ProyectoConfig, args) -> int:
    registro = Registro(cfg.raiz)
    r = registro.buscar(args.id)
    sched = obtener_scheduler(r.scheduler)
    try:
        print(sched.cancel(r.job_id))
    except SchedulerError as exc:
        print(f"No se pudo cancelar: {exc}")
        return 1
    registro.actualizar_estado(r.id, "CANCELADO")
    return 0


# ---------------------------------------------------------------------------
#  comparar
# ---------------------------------------------------------------------------
def cmd_comparar(cfg: ProyectoConfig, args) -> int:
    compuestos = cfg.resolver_compuestos(args.compuesto)
    tiempos = cfg.resolver_tiempos(args.tiempo) if args.tiempo else None

    filas = []
    for compuesto in compuestos:
        etiquetas = ([etiqueta_tiempo(t) for t in tiempos] if tiempos
                     else _tiempos_existentes(cfg, compuesto))
        for etq in etiquetas:
            ns = float(etq[:-2]) if etq.endswith("ns") else float(etq)
            layout = _layout(cfg, compuesto, ns)
            manifiesto = _leer_json(layout.manifiesto)
            filas.append([
                compuesto.clave,
                etq,
                f"{layout.segmentos_hechos()}/"
                f"{int(cfg.parametros_tiempo(ns)['segmentos'])}",
                "sí" if (layout.mmgbsa / "FINAL_MMGBSA.dat").exists() else "no",
                "sí" if (layout.ml / "features.csv").exists() else "no",
                (manifiesto or {}).get("huella", "-"),
                ((manifiesto or {}).get("entorno") or {}).get("servidor", "-"),
            ])
    if not filas:
        print("No hay ejecuciones que comparar todavía.")
        return 0
    print("Ejecuciones por compuesto y tiempo de simulación:\n")
    print(_tabla(filas, ["compuesto", "tiempo", "producción", "mmgbsa", "ml",
                         "huella", "servidor"]))
    print("\n  La huella resume los parámetros científicos: misma huella = "
          "comparable\n  entre servidores. Distinto tiempo de simulación cambia "
          "la huella a propósito.")
    return 0


def _reporte_tex(cfg: ProyectoConfig, args, compuestos) -> int:
    """Escribe las secciones generadas dentro del documento del paper."""
    from .reporte_tex import (bloque_comparacion, bloque_datos, cargar_corrida,
                              escribir_en_documento)

    documento = Path(args.documento)
    if not documento.is_absolute():
        documento = cfg.raiz / documento
    if not documento.exists():
        print(f"No existe el documento: {documento}")
        print("Indica otro con --documento.")
        return 1

    corridas = []
    for compuesto in compuestos:
        etiquetas = ([etiqueta_tiempo(t) for t in cfg.resolver_tiempos(args.tiempo)]
                     if args.tiempo else _tiempos_existentes(cfg, compuesto))
        for etq in etiquetas:
            limpio = etq.replace(" (heredado)", "")
            try:
                ns = float(limpio[:-2]) if limpio.endswith("ns") else float(limpio)
            except ValueError:
                continue
            corrida = cargar_corrida(compuesto, _layout(cfg, compuesto, ns))
            if corrida:
                corridas.append(corrida)

    if not corridas:
        print("No hay ninguna corrida con análisis (falta resumen.json).")
        print("Ejecuta antes la etapa 'analisis' del pipeline de resultados.")
        return 1

    for c in corridas:
        ml = "con rescoring ML" if c.ml else "sin rescoring ML"
        print(f"  {c.nombre} @ {c.etiqueta}: {ml}")

    bloques = {
        "datos_50ns": bloque_datos(corridas),
        "comparacion_md_ml": bloque_comparacion(corridas),
    }
    acciones = escribir_en_documento(documento, bloques,
                                     copia_seguridad=not args.sin_respaldo)
    respaldo = acciones.pop("_respaldo", None)
    print(f"\n{documento}")
    for nombre, accion in acciones.items():
        print(f"  bloque '{nombre}': {accion}")
    if respaldo:
        print(f"  copia previa: {respaldo}")
    print("\nSolo se reemplaza lo que hay entre las marcas % <<<dmlrun:...>>>;")
    print("el resto del documento queda intacto. Vuelve a compilar con:")
    print(f"    cd {documento.parent} && latexmk -pdf {documento.name}")
    return 0


def cmd_reporte(cfg: ProyectoConfig, args) -> int:
    from .reporte import escribir_indice, escribir_reporte, recolectar

    compuestos = cfg.resolver_compuestos(args.compuesto)

    if getattr(args, "formato", "html") in ("tex", "ambos"):
        codigo = _reporte_tex(cfg, args, compuestos)
        if args.formato == "tex":
            return codigo
        print()

    # La carpeta global de figuras (Img/) ha vivido en distintos niveles segun
    # como se clonara el proyecto; se usa la primera que exista.
    raiz_img = next((c for c in (cfg.raiz / "Img",
                                 cfg.raiz.parent / "Img",
                                 cfg.raiz.parent.parent / "Img") if c.is_dir()),
                    None)
    destino_base = (Path(args.destino) if args.destino
                    else cfg.raiz / "reportes")
    if not destino_base.is_absolute():
        destino_base = cfg.raiz / destino_base

    hechos: List = []
    sin_datos: List[str] = []

    for compuesto in compuestos:
        etiquetas = ([etiqueta_tiempo(t) for t in cfg.resolver_tiempos(args.tiempo)]
                     if args.tiempo else _tiempos_existentes(cfg, compuesto))
        for etq in etiquetas:
            limpio = etq.replace(" (heredado)", "")
            try:
                ns = float(limpio[:-2]) if limpio.endswith("ns") else float(limpio)
            except ValueError:
                continue
            layout = _layout(cfg, compuesto, ns)
            rec = recolectar(compuesto, layout, raiz_img=raiz_img)
            if rec.vacia:
                sin_datos.append(f"{compuesto.clave} @ {limpio}")
                continue
            carpeta = destino_base / f"{compuesto.clave}_{limpio}"
            info = escribir_reporte(rec, carpeta)
            hechos.append((rec, carpeta))
            print(f"[{compuesto.clave} · {limpio}] {info['figuras']} figura(s), "
                  f"{info['datos']} tabla(s)")
            print(f"  {info['informe']}")

    if not hechos:
        print("No encontré figuras ni resúmenes que recolectar.")
        if sin_datos:
            print("  sin datos: " + ", ".join(sin_datos))
        print("\nEjecuta primero las etapas de análisis, por ejemplo:")
        print("    dmlrun run --pipeline resultados --compuesto todos --tiempo 50 \\")
        print("           --etapas rmsd,analisis,mmgbsa --servidor hpc-cpu")
        return 1

    indice = escribir_indice(hechos, destino_base / "indice.html")
    print(f"\nÍndice con la comparación: {indice}")
    if sin_datos:
        print("  (sin datos todavía: " + ", ".join(sin_datos) + ")")
    print("\nLa carpeta es autocontenida: cada informe.html lleva las imágenes")
    print("incrustadas, así que se puede copiar, comprimir o enviar por correo.")
    return 0


def _tiempos_existentes(cfg, compuesto) -> List[str]:
    etiquetas = []
    if compuesto.runs_dir.exists():
        etiquetas += sorted(d.name for d in compuesto.runs_dir.iterdir() if d.is_dir())
    for d in sorted(compuesto.base_resultados.glob("md_*")):
        etq = "0.1ns" if d.name == "md_prueba" else d.name[3:]
        if etq not in etiquetas:
            etiquetas.append(etq)
    # Carpetas planas: resultados/10nsMacaridine, resultados/20nsMacaridine...
    sufijo = _nombre_entrega(compuesto, "")
    patron = re.compile(rf"^([\d.]+ns){re.escape(sufijo)}(_\d{{8}}-\d{{4}})?$")
    for d in sorted(_carpetas_de_compuesto(cfg), key=lambda x: x.name):
        m = patron.match(d.name)
        if m and m.group(1) not in etiquetas:
            etiquetas.append(m.group(1))
    return etiquetas


def _leer_json(ruta: Path) -> Optional[Dict]:
    try:
        return json.loads(Path(ruta).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


# ---------------------------------------------------------------------------
#  Parser
# ---------------------------------------------------------------------------
def construir_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="dmlrun",
        description="Orquestación de los pipelines de dockingML sobre HPC "
                    "(Slurm) y otros servidores Linux.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Ejemplos:
  dmlrun servidores
  dmlrun config --servidor hpc
  dmlrun run --pipeline resultados --compuesto todos --tiempo 10,20,50 --servidor hpc
  dmlrun run --pipeline ml --compuesto pallidine --tiempo 50 --fuente md --servidor hpc-cpu
  dmlrun run --pipeline ambos --compuesto todos --tiempo 20 --servidor linux-gpu
  dmlrun status
  dmlrun logs --id pallidine-resultados-50ns-20261105-094500 -f
""")
    p.add_argument("--raiz", default=None,
                   help="raíz del repositorio (por defecto se deduce; "
                        "también sirve $DMLRUN_RAIZ)")
    sub = p.add_subparsers(dest="comando", required=True)

    sub.add_parser("servidores", help="lista los perfiles de servidor")
    sub.add_parser("compuestos", help="lista los compuestos y sus ejecuciones")

    c = sub.add_parser("config", help="muestra y comprueba un perfil de servidor")
    c.add_argument("--servidor", default="auto")
    c.add_argument("--json", action="store_true")

    r = sub.add_parser("run", help="envía ejecuciones",
                       formatter_class=argparse.RawDescriptionHelpFormatter)
    r.add_argument("--pipeline", required=True,
                   help="resultados | ml | ambos")
    r.add_argument("--compuesto", required=True,
                   help="clave, lista separada por comas, o 'todos'")
    r.add_argument("--tiempo", default="50",
                   help="ns: '10', '10,20,50', 'prueba' o 'todos' (por defecto 50)")
    r.add_argument("--servidor", default="auto",
                   help="perfil de config/servidores.yaml; 'auto' lo detecta")
    r.add_argument("--etapas", default=None,
                   help="subconjunto de etapas separadas por comas "
                        "(resultados: sistema,min,cal,eq,prod,rmsd,analisis,mmgbsa; "
                        "ml: centro,docking,complejos,features,rescore)")
    r.add_argument("--layout", default="auto", choices=["auto", "runs", "legacy"],
                   help="auto reutiliza carpetas antiguas si existen")
    r.add_argument("--nueva-corrida", action="store_true",
                   help="estrena carpeta fechada en vez de reanudar la última")
    r.add_argument("--forzar", action="store_true",
                   help="repite etapas ya completadas")
    r.add_argument("--dry-run", action="store_true",
                   help="muestra el plan sin enviar nada")
    r.add_argument("--mostrar-script", action="store_true",
                   help="con --dry-run, imprime el script completo")
    r.add_argument("--solo-script", action="store_true",
                   help="escribe el script en estado/ sin enviarlo")
    r.add_argument("--ignorar-requisitos", action="store_true",
                   help="envía aunque falten entradas previas")
    # recursos (sobrescriben el perfil)
    r.add_argument("--cpus", type=int, default=None)
    r.add_argument("--memoria", type=int, default=None, help="GB")
    r.add_argument("--horas", type=int, default=None)
    # parámetros de MD
    r.add_argument("--t-calentamiento", type=float, default=None, help="ps")
    r.add_argument("--t-equilibrado", type=float, default=None, help="ps")
    r.add_argument("--segmentos", type=int, default=None)
    r.add_argument("--guardar-cada", type=float, default=None, help="ps por fotograma")
    r.add_argument("--max-segmentos", type=int, default=None,
                   help="correr como máximo N segmentos en este envío")
    r.add_argument("--frames-mmgbsa", type=int, default=None)
    r.add_argument("--semilla", type=int, default=2022,
                   help="semilla del integrador (reproducibilidad)")
    # parámetros de ML
    r.add_argument("--fuente", default="docking", choices=["docking", "md"],
                   help="poses para el pipeline de ML")
    r.add_argument("--motor", default="vina",
                   choices=["vina", "smina", "gnina", "qvina", "vina-cuda"])
    r.add_argument("--modelo", default=None, help="modelo .pkl para el rescoring")
    r.add_argument("--frames-ml", type=int, default=20,
                   help="fotogramas de MD a usar con --fuente md")

    s = sub.add_parser("status", help="estado de las ejecuciones")
    s.add_argument("--id", default=None, help="id completo o prefijo")
    s.add_argument("--compuesto", default=None)
    s.add_argument("--pipeline", default=None)
    s.add_argument("-n", type=int, default=15, help="últimas N ejecuciones")
    s.add_argument("--todas", action="store_true")

    l = sub.add_parser("logs", help="muestra el log de una ejecución")
    l.add_argument("--id", required=True)
    l.add_argument("-n", type=int, default=60)
    l.add_argument("-f", "--seguir", action="store_true")

    k = sub.add_parser("cancel", help="cancela una ejecución")
    k.add_argument("--id", required=True)

    comp = sub.add_parser("comparar",
                          help="compara tiempos de simulación de un compuesto")
    comp.add_argument("--compuesto", default="todos")
    comp.add_argument("--tiempo", default=None)

    rep = sub.add_parser(
        "reporte",
        help="reúne figuras y datos en un informe HTML autocontenido",
        description="Recoge las figuras y las tablas de cada corrida y escribe "
                    "un informe.html con las imágenes incrustadas, su "
                    "interpretación y los CSV al lado, listo para pasárselo a "
                    "alguien que no lanzó la simulación.")
    rep.add_argument("--compuesto", default="todos")
    rep.add_argument("--tiempo", default=None,
                     help="ns o lista; por defecto, todas las corridas que existan")
    rep.add_argument("--destino", default=None,
                     help="carpeta de salida (por defecto reportes/ en la raíz)")
    rep.add_argument("--formato", default="html", choices=["html", "tex", "ambos"],
                     help="html: informe autocontenido; tex: escribe las secciones "
                          "generadas dentro del documento del paper")
    rep.add_argument("--documento", default="paper/simulacion.tex",
                     help="documento LaTeX donde insertar las secciones (--formato tex)")
    rep.add_argument("--sin-respaldo", action="store_true",
                     help="no guardar una copia .bak del documento antes de escribir")

    return p


COMANDOS = {
    "servidores": cmd_servidores,
    "compuestos": cmd_compuestos,
    "config": cmd_config,
    "run": cmd_run,
    "status": cmd_status,
    "logs": cmd_logs,
    "cancel": cmd_cancel,
    "comparar": cmd_comparar,
    "reporte": cmd_reporte,
}


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = construir_parser().parse_args(argv)
    try:
        cfg = cargar_proyecto(args.raiz)
    except ConfigError as exc:
        print(f"Error de configuración: {exc}", file=sys.stderr)
        return 2
    try:
        return COMANDOS[args.comando](cfg, args)
    except (ConfigError, SchedulerError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    except KeyError as exc:
        print(f"Error: {exc.args[0] if exc.args else exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:  # pragma: no cover
        print("\nInterrumpido.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
