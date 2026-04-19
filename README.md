# Numerical Fluid Dynamics: Viscoplastic Flow & Stochastic Solvers

Professional implementation of numerical methods for complex fluid mechanics and physics problems, developed during research at the **Mechanical and Mathematical Faculty of MSU (Lomonosov Moscow State University)**.

## 🚀 Key Projects

### 1. Viscoplastic Flow Solvers (Bingham Fluid)
- **Problem:** Solving flow in channels and 2D flat domains (including L-shaped and cavity domains) using regularization and optimization methods.
- **Algorithms:** Implementation of **ALG2/ALG** iterative algorithms coupled with the **FISTA** (Fast Iterative Shrinkage-Thresholding Algorithm) for accelerated convergence.
- **Technology:** Developed in **C++** focusing on high-performance computations, sparse matrix handling, and Conjugate Gradient methods.
- **Numerical methods:** 3-point and 5-point sweep methods, iterative PNM/PTM schemes.

### 2. Atomic Interaction Modeling
- **Problem:** Stochastic simulation of hydrogen atoms and protons interactions.
- **Methodology:** Implemented **Monte Carlo** methods on C++ for physical process modeling.
- **Analysis:** Post-processing and data visualization using Python.

## 🛠 Tech Stack
- **Languages:** C++, Python, SQL, LaTeX
- **Academic Focus:** Numerical Methods, Optimization, Aeromechanics, Mathematical Statistics

## 📁 Repository Structure
- `Solver_L_ALG_FISTA/` — Flat L-shaped domain solver (FISTA, Sparse matrices).
- `Solver_ALG_FISTA_caverna/` — Cavity flow solver (FISTA, Sparse matrices).
- `solve_ALG2/` — Channel flow solvers (ALG2, L-shaped and square).
- `Monte_Karlo/` — Stochastic modeling of atomic interactions.

---
*Created by Ivan [Surname] as part of coursework at the Chair of Aeromechanics, MSU.*
