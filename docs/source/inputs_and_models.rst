Inputs and Models
=================

This page describes the model components that define a GPR forward simulation.

Model domain
------------

The domain is defined with physical extents and optional air and earth layers.
The Python I/O layer supports:

- air-only domains
- layered earth models
- rectangular and spherical anomalies
- receiver geometries in inline, endfire, and oblique directions

Layers and materials
--------------------

Each earth layer is defined by thickness, permittivity, conductivity, relative permeability, and magnetic conductivity.
The air layer is always present and is assigned a default dielectric permittivity and conductivity.

Anomalies
---------

Anomalies can be created as:

- `BoxAnomaly`
- `SphereAnomaly`

They are configured with spatial coordinates and electromagnetic properties.

Note: The old flat-parameter anomaly interface (e.g., `anomaly_x`, `anomaly_y`, `anomaly_z`,
and `anomaly_properties`) is still supported for backward compatibility but is considered
legacy. Prefer supplying a list of `BoxAnomaly`/`SphereAnomaly` objects to `GPRSurvey.build()`.

Source and receivers
--------------------

The source is configured by:

- antenna position
- current direction
- source type and waveform frequency content
- source moment and number of segments

Source refinement box: the optional source refinement box is implemented but marked as
experimental (see source code TODO in `io/inputs/sources.py`). It may be removed or
reworked in future versions; use with caution.

Receivers are placed in the domain using the receiver array definitions. The I/O layer supports inline, endfire, and oblique receiver placements.

Solver and PML
--------------

The solver configuration includes:

- solver type
- refinement settings
- maximum unknowns
- error tolerance
- output fields and VTK settings

Perfectly matched layers (PMLs) are configured in the I/O layer to absorb outgoing waves at the domain boundaries.
