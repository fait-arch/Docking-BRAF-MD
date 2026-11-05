# -*- coding: utf-8 -*-
"""
Los dos pipelines, descritos como etapas independientes del entorno.

* :class:`PipelineResults` — preparación del sistema, dinámica molecular
  (minimización → calentamiento → equilibrado → producción), RMSD, análisis y
  MM/GBSA. Es el camino que produce los resultados del artículo.
* :class:`MachineLearningPipeline` — docking, cálculo de descriptores de
  interacción y rescoring con el modelo de ML (``dockml.modern``).

Las dos clases comparten tres propiedades que son el objetivo de este módulo:

1. **Mismos compuestos**: el compuesto entra como dato (``Compuesto``), no
   codificado en el script. Un compuesto nuevo es una entrada en
   ``config/compuestos.yaml``.
2. **Cualquier tiempo de simulación**: el tiempo entra como ``RunLayout``, así
   que 10, 20 y 50 ns conviven en el mismo compuesto sin pisarse, y el pipeline
   de ML puede leer la ventana de MD que se le indique.
3. **Cualquier servidor**: los comandos no dependen del gestor de colas; el
   ``PerfilServidor`` sólo aporta rutas, plataforma de OpenMM y paralelismo.
"""
from __future__ import annotations

import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .config import Compuesto, PerfilServidor, ProyectoConfig
from .layout import RunLayout


@dataclass
class Etapa:
    """Una etapa ejecutable del pipeline."""

    nombre: str
    descripcion: str
    comandos: List[str] = field(default_factory=list)
    #: si es True, se omite cuando su salida ya existe (reanudable)
    omitible: bool = True
    #: razón por la que se omitió (la rellena ``Pipeline.plan``)
    omitida_por: Optional[str] = None

    @property
    def omitida(self) -> bool:
        return self.omitida_por is not None


class Pipeline:
    """Base común: resuelve rutas y compone la lista de etapas."""

    #: nombre con el que se selecciona desde la CLI
    clave: str = "base"
    #: etapas en orden
    ETAPAS: Sequence[str] = ()

    def __init__(self,
                 proyecto: ProyectoConfig,
                 compuesto: Compuesto,
                 layout: RunLayout,
                 servidor: PerfilServidor,
                 opciones: Optional[Dict] = None):
        self.proyecto = proyecto
        self.compuesto = compuesto
        self.layout = layout
        self.servidor = servidor
        self.opciones = dict(opciones or {})
        self.scripts = proyecto.raiz / "simulaciones" / "scripts"

    # -- utilidades ---------------------------------------------------------
    @property
    def python(self) -> str:
        return "python"

    def _script(self, nombre: str) -> str:
        """Ruta al script compartido, con reserva a la copia del compuesto."""
        compartido = self.scripts / nombre
        if compartido.exists():
            return shlex.quote(str(compartido))
        propio = self.compuesto.directorio / "scripts" / nombre
        return shlex.quote(str(propio))

    @staticmethod
    def _q(valor) -> str:
        return shlex.quote(str(valor))

    def etiqueta_trabajo(self) -> str:
        return f"{self.compuesto.clave}-{self.clave}-{self.layout.etiqueta}"

    # -- plan ---------------------------------------------------------------
    def etapas(self) -> List[Etapa]:  # pragma: no cover
        raise NotImplementedError

    def plan(self, seleccion: Optional[Sequence[str]] = None,
             forzar: bool = False) -> List[Etapa]:
        """Etapas a ejecutar, marcando las ya completadas como omitidas."""
        todas = {e.nombre: e for e in self.etapas()}
        desconocidas = set(seleccion or ()) - set(todas)
        if desconocidas:
            raise ValueError(
                f"etapas desconocidas para el pipeline '{self.clave}': "
                f"{', '.join(sorted(desconocidas))}. "
                f"Válidas: {', '.join(todas)}"
            )
        nombres = list(seleccion) if seleccion else [e.nombre for e in self.etapas()]
        plan: List[Etapa] = []
        for n in nombres:
            etapa = todas[n]
            if not forzar and etapa.omitible:
                razon = self._ya_hecha(n)
                if razon:
                    etapa.omitida_por = razon
            plan.append(etapa)
        return plan

    def _ya_hecha(self, etapa: str) -> Optional[str]:  # pragma: no cover
        return None

    def comandos(self, seleccion: Optional[Sequence[str]] = None,
                 forzar: bool = False) -> List[str]:
        """Las líneas de shell del trabajo, con cabeceras por etapa."""
        out: List[str] = [
            f'echo "== compuesto : {self.compuesto.nombre} ({self.compuesto.pubchem})"',
            f'echo "== pipeline  : {self.clave}"',
            f'echo "== tiempo    : {self.layout.etiqueta}"',
            f'echo "== servidor  : {self.servidor.nombre} ({self.servidor.scheduler})"',
            "",
        ]
        ejecutadas = 0
        for etapa in self.plan(seleccion, forzar=forzar):
            if etapa.omitida:
                out.append(f'echo "-- [{etapa.nombre}] omitida: {etapa.omitida_por}"')
                continue
            out.append(f'echo ""')
            out.append(f'echo "-- [{etapa.nombre}] {etapa.descripcion}"')
            out += etapa.comandos
            ejecutadas += 1
        if ejecutadas == 0:
            out.append('echo "Nada que hacer: todas las etapas estaban completas '
                       '(usa --forzar para repetirlas)."')
        return out


# ---------------------------------------------------------------------------
#  Pipeline de resultados: MD + análisis + MM/GBSA
# ---------------------------------------------------------------------------
class PipelineResults(Pipeline):
    """Dinámica molecular y análisis, parametrizado por tiempo de simulación.

    Etapas: ``sistema``, ``min``, ``cal``, ``eq``, ``prod``, ``rmsd``,
    ``analisis``, ``mmgbsa``. Las de dinámica son reanudables: si el segmento de
    producción ya tiene su ``.xml``, se salta.
    """

    clave = "resultados"
    ETAPAS = ("sistema", "min", "cal", "eq", "prod", "rmsd", "analisis", "mmgbsa")

    # -- parámetros derivados del tiempo ------------------------------------
    @property
    def params(self) -> Dict:
        p = self.proyecto.parametros_tiempo(self.layout.ns)
        # las opciones de la CLI mandan sobre los valores derivados
        for clave in ("t_calentamiento_ps", "t_equilibrado_ps", "segmentos",
                      "guardar_cada_ps", "frames_mmgbsa"):
            if self.opciones.get(clave) is not None:
                p[clave] = self.opciones[clave]
        return p

    # -- etapas -------------------------------------------------------------
    def etapas(self) -> List[Etapa]:
        L, C, S = self.layout, self.compuesto, self.servidor
        p = self.params
        q = self._q

        receptor = C.receptor_path(self.proyecto.raiz)
        ligando = C.ligando_path(self.proyecto.raiz)

        cmd_sistema = [
            f"{self.python} {self._script('01_preparar_sistema.py')} \\",
            f"    --receptor {q(receptor)} \\",
            f"    --ligando {q(ligando)} \\",
            f"    --out {q(L.sistema)} \\",
            f"    --residuos {q(C.receptor_residuos)} \\",
            f"    --ph {C.ph} \\",
            f"    --carga-ligando {int(C.carga_ligando)} \\",
            f"    --exhaustiveness {int(C.exhaustiveness)} \\",
            f"    --caja {C.caja_solvatacion} \\",
            f"    --sal {C.sal_molar}",
        ]

        def md(etapas_md: str) -> List[str]:
            flags = [
                f"--sistema {q(L.sistema)}",
                f"--out {q(L.md)}",
                f"--etapas {q(etapas_md)}",
                f"--t-calentamiento {p['t_calentamiento_ps']}",
                f"--t-equilibrado {p['t_equilibrado_ps']}",
                f"--t-produccion {p['t_produccion_ps']}",
                f"--segmentos {int(p['segmentos'])}",
                f"--guardar-cada {p['guardar_cada_ps']}",
                f"--plataforma {q(S.plataforma_openmm)}",
                f"--hilos {int(S.cpus)}",
                f"--semilla {int(self.opciones.get('semilla', 2022))}",
            ]
            if self.opciones.get("max_segmentos"):
                flags.append(f"--max-segmentos {int(self.opciones['max_segmentos'])}")
            lineas = [f"{self.python} {self._script('02_dinamica.py')} \\"]
            lineas += [f"    {f} \\" for f in flags]
            lineas.append(f'    2>&1 | tee -a {q(L.md / "run.log")}')
            return lineas

        cmd_rmsd = [
            f"{self.python} {self._script('03_rmsd.py')} \\",
            f"    --sistema {q(L.sistema)} \\",
            f"    --md {q(L.md)} \\",
            f"    --out {q(L.resultados)} \\",
            f"    --titulo {q(C.titulo)}",
        ]

        cmd_analisis = [
            f"{self.python} {self._script('04_analisis_reporte.py')} \\",
            f"    --base {q(C.directorio)} \\",
            f"    --md {q(L.rel(L.md))} \\",
            f"    --fig {q(L.resultados)} \\",
            f"    --out {q(L.analisis)}",
        ]

        return [
            Etapa("sistema",
                  "Preparación del sistema (docking Vina + topología AMBER)",
                  cmd_sistema),
            Etapa("min", "Minimización (10 000 it., k = 10 kcal/mol/Å²)", md("min")),
            Etapa("cal", f"Calentamiento 0 → 310 K en {p['t_calentamiento_ps']} ps (NVT)",
                  md("cal")),
            Etapa("eq", f"Equilibrado {p['t_equilibrado_ps']} ps (NPT, 310 K, 1 atm)",
                  md("eq")),
            Etapa("prod",
                  f"Producción {self.layout.ns:g} ns en {int(p['segmentos'])} "
                  f"segmento(s) reanudable(s)",
                  md("prod")),
            Etapa("rmsd", "RMSD de Cα y del ligando", cmd_rmsd, omitible=False),
            Etapa("analisis", "Análisis completo y figuras", cmd_analisis, omitible=False),
            Etapa("mmgbsa", f"MM/GBSA ({int(p['frames_mmgbsa'])} fotogramas)",
                  self._cmd_mmgbsa(), omitible=True),
        ]

    def _cmd_mmgbsa(self) -> List[str]:
        """MM/GBSA sobre la trayectoria del tiempo elegido.

        El intervalo de fotogramas se calcula en tiempo de ejecución a partir de
        los .dcd que existan, no del tiempo nominal: así funciona igual si la
        producción se cortó antes de terminar.
        """
        L, C = self.layout, self.compuesto
        p = self.params
        q = self._q
        return [
            f"mkdir -p {q(L.mmgbsa)}",
            f"{self.python} {self._script('05_preparar_mmgbsa.py')} \\",
            f"    --sistema {q(L.sistema)} \\",
            f"    --md {q(L.md)} \\",
            f"    --out {q(L.mmgbsa)} \\",
            f"    --frames {int(p['frames_mmgbsa'])} \\",
            f"    --guardar-cada {p['guardar_cada_ps']} \\",
            f"    --titulo {q(f'MM/GBSA {C.nombre}-BRAF V600E ({L.etiqueta})')}",
            "",
            f"DCDS=$(ls {q(str(L.md))}/04_prod_*.dcd | sort | tr '\\n' ' ')",
            f'if [ -z "$DCDS" ]; then echo "sin trayectorias de producción en '
            f'{L.md}"; exit 1; fi',
            f"( cd {q(L.mmgbsa)} && MMPBSA.py -O -i mmgbsa.in \\",
            "      -o FINAL_MMGBSA.dat -do FINAL_DECOMP.dat \\",
            "      -cp com.prmtop -rp rec.prmtop -lp lig.prmtop -y $DCDS \\",
            "      -eo energias_por_frame.csv -deo decomp_por_frame.csv \\",
            "      > mmpbsa.log 2>&1 )",
            "",
            f"{self.python} {self._script('05_mmgbsa_figura.py')} \\",
            f"    --mmgbsa {q(L.mmgbsa)} \\",
            f"    --fig {q(L.resultados)} \\",
            f'    --dt-ps "$(cat {q(L.mmgbsa / "dt_ps.txt")})"',
        ]

    # -- reanudación --------------------------------------------------------
    def _ya_hecha(self, etapa: str) -> Optional[str]:
        L = self.layout
        if etapa == "sistema" and self.compuesto.sistema_preparado:
            return "sistema/sistema.prmtop ya existe"
        hechas = L.etapas_completadas()
        if etapa in ("min", "cal", "eq") and hechas.get(etapa):
            return f"{L.rel(L.md)} ya tiene el estado de '{etapa}'"
        if etapa == "prod":
            total = int(self.params["segmentos"])
            hechos = L.segmentos_hechos()
            if hechos >= total:
                return f"los {total} segmentos de producción ya están hechos"
        if etapa == "mmgbsa" and (L.mmgbsa / "FINAL_MMGBSA.dat").exists():
            return "FINAL_MMGBSA.dat ya existe"
        return None


# ---------------------------------------------------------------------------
#  Pipeline de machine learning
# ---------------------------------------------------------------------------
class MachineLearningPipeline(Pipeline):
    """Docking, descriptores de interacción y rescoring con ML.

    Dos fuentes de poses, seleccionables con ``--fuente``:

    * ``docking`` — redocking con el motor elegido sobre el receptor preparado.
      No depende del tiempo de simulación; es la referencia común.
    * ``md`` — fotogramas representativos de la producción del tiempo indicado
      (10, 20 o 50 ns). Permite comparar cómo cambia el rescoring según la
      ventana simulada, sobre los mismos compuestos.
    """

    clave = "ml"
    ETAPAS = ("centro", "docking", "complejos", "features", "rescore")

    @property
    def fuente(self) -> str:
        return str(self.opciones.get("fuente", "docking"))

    @property
    def motor(self) -> str:
        return str(self.opciones.get("motor", "vina"))

    @property
    def modelo(self) -> Optional[str]:
        return self.opciones.get("modelo")

    def etapas(self) -> List[Etapa]:
        L, C, S = self.layout, self.compuesto, self.servidor
        q = self._q
        n_frames = int(self.opciones.get("frames_ml", 20))

        def entrada(paso: str) -> List[str]:
            flags = [
                f"--paso {paso}",
                f"--fuente {q(self.fuente)}",
                f"--sistema {q(L.sistema)}",
                f"--md {q(L.md)}",
                f"--out {q(L.ml)}",
                f"--frames {n_frames}",
                f"--ligando-resname {q(C.ligando_resname)}",
                f"--bolsillo {q(','.join(str(r) for r in C.bolsillo_residuos))}",
            ]
            lineas = [f"{self.python} {self._script('06_ml_entrada.py')} \\"]
            lineas += [f"    {f} \\" for f in flags[:-1]]
            lineas.append(f"    {flags[-1]}")
            return lineas

        # -- docking ---------------------------------------------------------
        cmd_docking = [
            f'CENTRO=$(cat {q(L.ml / "centro_caja.txt")})',
            "dml-dock \\",
            f"    --receptor {q(L.sistema / 'receptor.pdbqt')} \\",
            f"    --engine {q(self.motor)} \\",
            "    --center $CENTRO \\",
            f"    --size {' '.join(str(x) for x in C.caja_vina)} \\",
            f"    --exhaustiveness {int(C.exhaustiveness)} \\",
            f"    --ncpu {int(S.cpus)} \\",
            f"    --outdir {q(L.ml / 'docked')} \\",
            f"    {q(L.sistema / 'ligando.pdbqt')}",
        ]

        # -- descriptores ----------------------------------------------------
        cmd_features = [
            "dml-features \\",
            f"    --input {q(L.ml / 'complejos.txt')} \\",
            f"    --output {q(L.ml / 'features.csv')} \\",
            "    --onion --plif",
        ]

        # -- rescoring -------------------------------------------------------
        if self.modelo:
            cmd_rescore = [
                "dml-rescore \\",
                f"    --features {q(L.ml / 'features.csv')} \\",
                f"    --model {q(self.modelo)} \\",
                f"    --output {q(L.ml / 'rescored.csv')}",
            ]
            desc_rescore = f"Rescoring con el modelo {Path(str(self.modelo)).name}"
        else:
            cmd_rescore = [
                f'echo "Sin --modelo: se dejan los descriptores en '
                f'{L.rel(L.ml / "features.csv")} sin reescalar."',
                f'echo "Para entrenar uno:  dml-train -X {L.rel(L.ml / "features.csv")}'
                ' -y etiquetas.csv -o modelo.pkl"',
            ]
            desc_rescore = "Rescoring (omitido: no se indicó --modelo)"

        if self.fuente == "md":
            desc_complejos = (f"{n_frames} fotogramas de la producción de "
                              f"{L.etiqueta} como complejos")
        else:
            desc_complejos = "Complejos receptor-ligando a partir de las poses acopladas"

        etapas = [
            Etapa("centro", "Centro de la caja de docking (bolsillo de ATP)",
                  entrada("centro"), omitible=False),
            Etapa("docking", f"Docking con {self.motor}", cmd_docking),
            Etapa("complejos", desc_complejos, entrada("complejos"), omitible=False),
            Etapa("features", "Descriptores de interacción (OnionNet + PLIF)",
                  cmd_features, omitible=False),
            Etapa("rescore", desc_rescore, cmd_rescore, omitible=False),
        ]
        if self.fuente == "md":
            # los complejos vienen de la trayectoria: no hay nada que acoplar
            etapas = [e for e in etapas if e.nombre != "docking"]
        return etapas

    def _ya_hecha(self, etapa: str) -> Optional[str]:
        L = self.layout
        if etapa == "docking" and (L.ml / "docked" / "docking_summary.csv").exists():
            return "docking_summary.csv ya existe"
        return None

    def requisitos(self) -> List[str]:
        """Comprobaciones previas que la CLI muestra como aviso."""
        faltan = []
        if self.fuente == "md" and not self.layout.trayectorias():
            faltan.append(
                f"no hay trayectorias en {self.layout.rel(self.layout.md)}; "
                f"ejecuta primero el pipeline 'resultados' para {self.layout.etiqueta}"
            )
        if self.fuente == "docking" and not self.compuesto.sistema_preparado:
            faltan.append(
                f"falta {self.compuesto.clave}/sistema/sistema.prmtop; "
                "ejecuta la etapa 'sistema' del pipeline 'resultados'"
            )
        return faltan


PIPELINES = {
    PipelineResults.clave: PipelineResults,
    MachineLearningPipeline.clave: MachineLearningPipeline,
}


def construir_pipeline(clave: str,
                       proyecto: ProyectoConfig,
                       compuesto: Compuesto,
                       layout: RunLayout,
                       servidor: PerfilServidor,
                       opciones: Optional[Dict] = None) -> Pipeline:
    """Fábrica: ``'resultados'`` o ``'ml'``."""
    if clave not in PIPELINES:
        raise ValueError(
            f"pipeline '{clave}' desconocido. Válidos: {', '.join(PIPELINES)}"
        )
    return PIPELINES[clave](proyecto, compuesto, layout, servidor, opciones)
