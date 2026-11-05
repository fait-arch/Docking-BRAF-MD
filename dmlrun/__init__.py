# -*- coding: utf-8 -*-
"""
dmlrun — orquestación de los pipelines de dockingML sobre distintos servidores.

El paquete separa tres cosas que antes estaban mezcladas en los scripts:

* **qué** se ejecuta  -> :mod:`dmlrun.pipelines` (``PipelineResults`` y
  ``MachineLearningPipeline``), que describen las etapas de forma independiente
  del compuesto y del tiempo de simulación;
* **dónde** se ejecuta -> :mod:`dmlrun.config` + :mod:`dmlrun.schedulers`
  (perfiles de servidor en ``config/servidores.yaml``, backends Slurm y local);
* **dónde se guarda** -> :mod:`dmlrun.layout` (``runs/<tiempo>/``), con registro
  y manifiesto de reproducibilidad en :mod:`dmlrun.registry`.

Uso típico desde la terminal::

    dmlrun run --pipeline ambos --compuesto todos --tiempo 10,20,50 --servidor hpc
    dmlrun status
"""
from __future__ import annotations

__version__ = "1.0.0"

from .config import (  # noqa: F401
    Compuesto,
    PerfilServidor,
    ProyectoConfig,
    cargar_proyecto,
    detectar_servidor,
    encontrar_raiz,
)
from .layout import RunLayout, etiqueta_tiempo  # noqa: F401
from .pipelines import (  # noqa: F401
    Etapa,
    MachineLearningPipeline,
    Pipeline,
    PipelineResults,
    construir_pipeline,
)
from .registry import Registro, RunRecord  # noqa: F401
from .schedulers import JobSpec, LocalScheduler, SlurmScheduler, obtener_scheduler  # noqa: F401

__all__ = [
    "Compuesto",
    "PerfilServidor",
    "ProyectoConfig",
    "cargar_proyecto",
    "detectar_servidor",
    "encontrar_raiz",
    "RunLayout",
    "etiqueta_tiempo",
    "Etapa",
    "Pipeline",
    "PipelineResults",
    "MachineLearningPipeline",
    "construir_pipeline",
    "Registro",
    "RunRecord",
    "JobSpec",
    "LocalScheduler",
    "SlurmScheduler",
    "obtener_scheduler",
    "__version__",
]
