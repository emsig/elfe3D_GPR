# Outputs and Postprocessing

This page explains the file outputs that `elfe3D_GPR` writes and the tools available for postprocessing.

## Solver inputs and outputs

The Python I/O layer writes the following files for each simulation:

- `elfe3D_input.txt` — main solver configuration
- `source.txt` — source segment and excitation definition
- `regionparameters.txt` — material properties per region
- `<model>.poly` — TetGen geometry file
- solver output files in `out_<experiment_name>`

## Field outputs

The solver writes electric and magnetic field results to the configured output files. These are typically stored in:

- `out_<experiment_name>/electric_fields*`
- `out_<experiment_name>/magnetic_fields*`

## Python postprocessing tools

The repository includes helper modules in `io/outputs` for reading and visualizing results, including:

- `io.outputs.fieldreader`
- `io.outputs.visualize`
- `io.outputs.postprocess`

These modules are intended to support notebook-based analysis of solver output.
