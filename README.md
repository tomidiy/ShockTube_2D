# 2D Compressible Shock Tube Simulation

This repository contains a Python implementation of a **2D compressible shock tube simulation** using the **HLLC Riemann solver** to solve the Euler equations for inviscid compressible flow.
The simulation models the evolution of density, velocity, and pressure in a 2D domain, initialized with a shock tube setup (discontinuous left/right states).

## Features
- Solves the 2D Euler equations using a finite volume method.
- Implements the HLLC (Harten-Lax-van Leer-Contact) Riemann solver for accurate flux computation.
- Supports a 2D shock tube problem with customizable parameters (grid size, CFL number, etc.).
- Outputs `.npz` files (density, velocities, pressure) and PNG contour plots in the `results/` directory.
- Modular design for easy extension (e.g., to 3D or other Riemann solvers).

## Requirements
- Python 3.8+
- NumPy (>=1.21.0)
- Matplotlib (>=3.4.0)

## Installation
1. **Clone the repository**:
   ```bash
   git clone https://github.com/[your-username]/shock_tube_2d.git
   cd shock_tube_2d
   ```

2. **Install dependencies**:
   ```bash
      pip install -r requirements.txt
   ```

## Usage
- Run the simulation using the main script:
```bash
    python src/main.py
```
- Outputs: `.npz` files in results/ containing density, x-velocity, y-velocity, and pressure fields at specified intervals.
  PNG contour plots in `results/` visualizing the fields.

- **Default Parameters**:
  - Grid: 200x200 cells
  - Domain: 1.0x1.0 units
  - End time: 0.2
  - CFL number: 0.8
  - Output interval: 0.01
  - Specific heat ratio (γ): 1.4 (air)
  - Shock tube states: Left (ρ=1, u=0, v=0, p=1), Right (ρ=0.125, u=0, v=0, p=0.1)
- Modify parameters in `src/main.py` to adjust grid size, CFL, or simulation duration.



## Project Structure
```
shock_tube_2d/
├── src/
│   ├── main.py       # Main driver script
│   ├── grid.py      # Grid setup and initial conditions
│   ├── solver.py    # HLLC solver and flux calculations
│   ├── output.py    # Data saving and visualization
├── results/         # Output directory for .npz and .png files
├── requirements.txt # Dependencies
├── .gitignore       # Git ignore file
└── LICENSE          # MIT License
```

Example OutputAfter running, check `results/` for: 
-Data: `fields_t0.010.npz`, `fields_t0.020.npz`, etc., with arrays for `ρ`, `u`, `v`, `p`.
-Plots: `figure_t0.010.png`, etc., showing 2x2 contour plots of density, x-velocity, y-velocity, and pressure.

##Contributing
Contributions are welcome! Please:
- Fork the repository.
- Create a feature branch (git checkout -b feature-name).
- Commit changes (git commit -m "Add feature").
- Push to the branch (git push origin feature-name).
- Open a pull request.

##License
This project is licensed under the MIT License. See the LICENSE file for details.

##Acknowledgments
- Inspired by computational fluid dynamics literature, particularly Toro's `Riemann Solvers and Numerical Methods for Fluid Dynamics`.
- Built with NumPy and Matplotlib for numerical and visualization support.





