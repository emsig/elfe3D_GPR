"""
writeinputfiles.py
------------------
Writes all elfe3D FEM input files (everything except the .poly file).

FEMInputWriter takes the assembled domain objects and writes:
  - elfe3D_input.txt  (main solver configuration file)
  - source.txt        (source segment coordinates)
  - regionparameters.txt (material properties per region attribute)

The .poly file is handled separately by PolyAssembler in writetetgenpoly.py.
"""

import os
import textwrap
from pathlib import Path
from datetime import datetime

from .domain import ModelDomain
from .geolayers import LayerStack
from .sources import SourceAntenna
from .receivers import ReceiverArray
from .solver import SolverConfig
from .pml import PMLConfig


_COORD_DIGITS = 6

def _r(x: float) -> float:
    """Round a coordinate to standard output precision (6 decimal places)."""
    return round(float(x), _COORD_DIGITS)


class FEMInputWriter:
    """
    Writes the elfe3D FEM input files from assembled domain objects.

    Parameters
    ----------
    domain    : ModelDomain
    layers    : LayerStack
    source    : SourceAntenna
    receivers : ReceiverArray
    solver    : SolverConfig
    pml       : PMLConfig
    output_dir : Path: all files written here
    io_config  : dict: file name overrides (optional, see defaults below)

    Default filenames
    -----------------
    base_input_file    : "elfe3D_input.txt"
    source_file        : "source.txt"
    region_params_file : "regionparameters.txt"
    """

    # Default file names — override via io_config if needed
    _DEFAULTS = {
        "base_input_file":    "elfe3D_input.txt",
        "source_file":        "source.txt",
        "region_params_file": "regionparameters.txt",
    }

    def __init__(
        self,
        domain: ModelDomain,
        layers: LayerStack,
        source: SourceAntenna,
        receivers: ReceiverArray,
        solver_cfg: SolverConfig,
        pml: PMLConfig,
        output_dir: Path | str,
        model_file_path: str = "",   # path written into input file for .poly/.1
        output_E_file: str = "",
        output_H_file: str = "",
        io_config: dict | None = None,
    ):
        self.domain      = domain
        self.layers      = layers
        self.source      = source
        self.receivers   = receivers
        self.solver_cfg  = solver_cfg
        self.pml         = pml
        self.output_dir  = Path(output_dir)

        # File paths written into elfe3D_input.txt (not the path of that file)
        self.model_file_path = model_file_path
        self.output_E_file   = output_E_file
        self.output_H_file   = output_H_file

        cfg = io_config or {}
        self.base_input_file    = self.output_dir / cfg.get("base_input_file",    self._DEFAULTS["base_input_file"])
        self.source_file        = self.output_dir / cfg.get("source_file",        self._DEFAULTS["source_file"])
        self.region_params_file = self.output_dir / cfg.get("region_params_file", self._DEFAULTS["region_params_file"])


    def write_all(self, regions: list) -> None:
        """
        Write all FEM input files.

        Parameters
        ----------
        regions : assembled region list from PolyAssembler.regions
                  (needed for regionparameters.txt)
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.write_base_input_file()
        self.write_source_file()
        self.write_region_params_file(regions)

    # ------------------------------------------------------------------
    # elfe3D_input.txt
    # ------------------------------------------------------------------

    def write_base_input_file(self) -> Path:
        """
        Write the main elfe3D_input.txt file.
        """
        pml_buff = self.pml.total_thickness
        sc = self.solver_cfg

        with open(self.base_input_file, 'w') as f:
            f.write(f"solver                  {sc.solver}\n")
            f.write("model_size\n")
            f.write(f"{_r(self.domain.x_min - pml_buff)}\t{_r(self.domain.y_min - pml_buff)}\t{_r(self.domain.z_min - pml_buff)}\n")
            f.write(f"{_r(self.domain.x_max + pml_buff)}\t{_r(self.domain.y_max + pml_buff)}\t{_r(self.domain.z_max + pml_buff)}\n")

            f.write(f"num_freq                {len(self.source.f_list)}\n")
            for freq in self.source.f_list:
                f.write(f"{freq}\n")

            f.write(f"num_rec                 {self.receivers.num_receivers}\n")
            for x, y, z in zip(self.receivers.x, self.receivers.y, self.receivers.z):
                f.write(f"{_r(x)} {_r(y)} {_r(z)}\n")

            f.write(f"output_E_file           {self.output_E_file}\n")
            f.write(f"output_H_file           {self.output_H_file}\n")
            f.write(f"source_type             {self.source.source_type}\n")
            f.write(f"{_r(self.source.x_extents[0])} {_r(self.source.y_extents[0])} {_r(self.source.z_extents[0])}\n")
            f.write(f"{_r(self.source.x_extents[1])} {_r(self.source.y_extents[1])} {_r(self.source.z_extents[1])}\n")
            f.write(f"current_direction       {self.source.current_direction}\n")
            f.write(f"source_moment           {self.source.source_moment}\n")
            f.write(f"PEC_present             {0}\n")
            f.write(f"num_PEC                 {0}\n")
            f.write(f"model_file_name         {self.model_file_path}\n")
            f.write(f"maxRefSteps             {sc.max_ref_steps}\n")
            f.write(f"maxUnknowns             {int(sc.max_unknowns)}\n")
            f.write(f"betaRef                 {sc.beta_ref}\n")
            f.write(f"accuracyTol             {sc.accuracy_tol:.6f}\n")
            f.write(f"vtkRef                  {sc.vtk}\n")
            f.write(f"errorEst_method         {sc.error_est_method}\n")
            f.write(f"refStrategy             {sc.ref_strategy}\n")
            f.write(f"output_fields_vtk       {sc.output_fields_vtk}\n")
            f.write(f"PML_present             {int(self.pml.present)}\n")
            f.write(f"PML_thickness           {self.pml.layer_thickness}\n")
            f.write(f"PML_decay_type          {sc.pml_decay_type}\n")

        return self.base_input_file

    # ------------------------------------------------------------------
    # source.txt
    # ------------------------------------------------------------------

    def write_source_file(self) -> Path:
        """
        Write the source segment coordinate file.
        """
        with open(self.source_file, 'w') as f:
            f.write(f"{2}\n")
            for i in range(2):
                f.write(
                    f"{_r(self.source.x_extents[i])} "
                    f"{_r(self.source.y_extents[i])} "
                    f"{_r(self.source.z_extents[i])}\n"
                )
        return self.source_file

    # ------------------------------------------------------------------
    # regionparameters.txt
    # ------------------------------------------------------------------

    def write_region_params_file(self, regions: list) -> Path:
        """
        Write the region electromagnetic parameters file.

        Parameters
        ----------
        regions : list of 10-element region lists from PolyAssembler
                  region[0] = attribute, region[7] = rho,
                  region[8] = mu_r, region[9] = eps_r
        """
        with open(self.region_params_file, 'w') as f:
            f.write("# eleattr\n")
            f.write(f"{len(regions)}\n")
            f.write("# eleattr rho mu_r epsilon_r\n")
            for region in regions:
                f.write(f"{region[0]} {region[7]} {region[8]} {region[9]}\n")
        return self.region_params_file

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        return f"FEMInputWriter(output_dir={self.output_dir})"
