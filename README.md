# Numerical Fluid Dynamics: Viscoplastic Flow Solvers

Professional implementation of numerical methods for specialized fluid mechanics problems, developed during research at the **Mechanical and Mathematical Faculty of MSU (Lomonosov Moscow State University)**.

## 🚀 Key Projects

### 1. Flat L-Shaped Domain Solver
- **Problem:** Modeling viscoplastic flow in complex L-shaped geometries using regularization techniques.
- **Algorithms:** Implementation of **FISTA** (Fast Iterative Shrinkage-Thresholding Algorithm) for accelerated minimization of the viscoplastic functional.
- **Technology:** **C++** with automated multithreading (`ThreadPool.h`) and sparse matrix optimizations.
- **Visualization:** Integrated Python scripts for modeling flow fields and rigid zones.

### 2. Cavity Flow Solver (Caverna)
- **Problem:** Simulating 2D viscoplastic flow in a square cavity.
- **Optimization:** Comparison between standard ALM (Augmented Lagrangian Method) and accelerated FISTA variants.
- **Numerical methods:** High-precision discrete solvers for Poisson equations and velocity field updates.

### 3. Channel Flow Solvers (ALG2)
- **Problem:** Classical viscoplastic flow in simplified channel geometries.
- **Methodology:** Robust **ALG2** iterative algorithm for reliable convergence in viscoplastic modeling.

## 🛠 Tech Stack
- **Languages:** C++, Python, SQL, LaTeX
- **Numerical libraries:** Custom sparse solvers, Conjugate Gradient, Poisson solvers.
- **Academic Focus:** Numerical Methods, Non-Newtonian Fluids, Aeromechanics.

## 📁 Repository Structure
- `Solver_L_ALG_FISTA/` — L-shaped domain solver (Main project).
- `Solver_ALG_FISTA_caverna/` — Cavity flow solver implementation.
- `solve_ALG2/` — Collection of channel flow solvers and prototypes.

---
*Created by Ivan Platonychev as part of coursework at the Chair of Aeromechanics, MSU.*
