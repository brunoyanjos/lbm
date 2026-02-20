# LBM

CUDA C++ LBM framework for stencil analysis, boundary condition studies, and performance evaluation.

---

## Overview

This repository contains a research-oriented CUDA C++ implementation of the Lattice Boltzmann Method (LBM). The project is designed as a structured framework for systematic numerical and architectural studies.

The main research directions include:

- Stencil design and comparison (D2Q9 implemented, D2V17/D3Q19/D3Q27 planned)
- Boundary condition treatments (regularized formulations)
- Solver framework strategies (e.g., double-array, AB variants, and related designs)
- Performance benchmarking in MLUPS
- Progressive extension toward higher-order and advanced formulations

The project is currently under active development and includes preliminary results.

---

## Current Features

- 2D D2Q9 stencil
- BGK collision operator
- Moment-based evaluation
- Regularized boundary condition treatment
- Square cavity geometry
- Explicit warmup phase before performance measurement
- GPU timing via CUDA events
- Wall-clock timing via `std::chrono`
- MLUPS reporting (GPU and wall time)
- VTK output
- Performance logging
- Terminal progress bar with ETA and partial MLUPS

---

## Planned Extensions

- D2V17 stencil
- D3Q19 and D3Q27 stencils
- Higher-order formulations
- Nonlinear system solvers
- Additional geometries
- Extended output fields (velocity, kinetic energy, etc.)
- Alternative memory layouts and solver frameworks

---

## Numerical Model

The solver currently implements:

- BGK collision model
- Moment-based evaluation
- Regularized boundary reconstruction

The regularized boundary treatment follows formulations described in:

Hegele Jr., L. A., Scagliarini, A., Sbragaglia, M., Mattila, K. K., Philippi, P. C., Puleri, D. F., Gounley, J., & Randles, A.  
High-Reynolds-number turbulent cavity flow using the lattice Boltzmann method.  
Physical Review E, 98, 043302 (2018).  
https://doi.org/10.1103/PhysRevE.98.043302

The framework is structured to support further extensions toward higher-order moment representations and alternative collision strategies.

---

## Performance Measurement

Performance is reported in MLUPS (Million Lattice Updates Per Second).

Benchmark procedure:

- Explicit warmup phase before timing
- GPU time measured using CUDA events
- Wall time measured using `std::chrono`
- MLUPS computed as:

MLUPS = (NX × NY × measured_steps) / gpu_time / 1e6

Both GPU-based and wall-based MLUPS are reported.

Partial MLUPS values are displayed during runtime through the progress interface.

---

## Project Structure

src/
├── app/ # Simulation control, CUDA configuration, progress UI
├── core/ # Geometry, types, utilities, configuration
├── io/ # Output routines (VTK, metadata, performance)
└── lbm/ # Core LBM implementation
├── boundary/
├── collision/
├── domain/
├── moment_evaluation/
├── population/
├── state/
└── stencils/

- `app/` handles execution flow and benchmarking.
- `lbm/` contains the numerical method implementation.
- `io/` manages simulation outputs.
- `core/` provides shared definitions and utilities.

---

## Build

### Requirements

- Linux environment (tested on Ubuntu via WSL)
- CUDA Toolkit ≥ 12
- C++17-compatible compiler
- NVIDIA GPU with CUDA support

### Compile (default configuration)

./compile.sh

### Example with custom configuration

STENCIL=D2Q9 REAL=double DEBUG=1 ./compile.sh

---

## Run

The executable is generated under:
build/<configuration>/sim

By default, `compile.sh` runs the simulation automatically.

Outputs are written to:
runs/<timestamp>/
├── vtk/
└── logs/

---

## Status

This is a research codebase under active development.

Issues are currently disabled.  
Comments and suggestions are welcome.

---

## License

MIT License
