#include "LShapeSolver2D.h"
#define NOMINMAX
#include <windows.h>
#include <iostream>

// ============================================================
// Constructor / Domain Helpers
// ============================================================

LShapeSolver2D::LShapeSolver2D(LShapeConfig config) : cfg(config) {
    nx_max = std::max(cfg.nx1, cfg.nx2);
    ny_max = std::max(cfg.ny1, cfg.ny2);
    // Create pool ONCE with 2 pre-allocated threads (u and v solved in parallel)
    pool = std::make_unique<ThreadPool>(2);
    prepare_fields();
    precompute_matrices();
}

bool LShapeSolver2D::is_in_L_shape_u(int j, int i) {
    if (j >= 1 && j <= cfg.ny1 && i >= 1 && i <= cfg.nx2 - 1) return true;
    if (j >= cfg.ny1 + 1 && j <= cfg.ny2 && i >= 1 && i <= cfg.nx1 - 1) return true;
    return false;
}

bool LShapeSolver2D::is_in_L_shape_v(int j, int i) {
    if (j >= 1 && j <= cfg.ny1 - 1 && i >= 1 && i <= cfg.nx2) return true;
    if (j == cfg.ny1 && i >= 1 && i <= cfg.nx1) return true;
    if (j >= cfg.ny1 + 1 && j <= cfg.ny2 - 1 && i >= 1 && i <= cfg.nx1) return true;
    return false;
}

bool LShapeSolver2D::is_in_L_shape_p(int j, int i) {
    if (j >= 1 && j <= cfg.ny1 && i >= 1 && i <= cfg.nx2) return true;
    if (j >= cfg.ny1 + 1 && j <= cfg.ny2 && i >= 1 && i <= cfg.nx1) return true;
    return false;
}

std::string LShapeSolver2D::out_path(const std::string& name) {
    return cfg.output_dir + "\\" + name;
}

// ============================================================
// Prepare fields
// ============================================================

void LShapeSolver2D::prepare_fields() {
    int ny = ny_max + 2, nx = nx_max + 2;
    u.assign(ny, std::vector<double>(nx, 0.0));
    v.assign(ny, std::vector<double>(nx, 0.0));
    p.assign(ny, std::vector<double>(nx, 0.0));
    u_tilde.assign(ny, std::vector<double>(nx, 0.0));
    v_tilde.assign(ny, std::vector<double>(nx, 0.0));
    fu.assign(ny, std::vector<double>(nx, 0.0));
    fv.assign(ny, std::vector<double>(nx, 0.0));
    Ap.assign(ny, std::vector<double>(nx, 0.0));
    r1.assign(ny, std::vector<double>(nx, 0.0));
    p_prev.assign(ny, std::vector<double>(nx, 0.0));
    F.assign(ny, std::vector<double>(nx, 0.0));

    // Bingham allocations
    f11.assign(ny, std::vector<double>(nx, 0.0));
    f12.assign(ny, std::vector<double>(nx, 0.0));
    f22.assign(ny, std::vector<double>(nx, 0.0));

    tau11.assign(ny, std::vector<double>(nx, 0.0));
    tau12.assign(ny, std::vector<double>(nx, 0.0));
    tau22.assign(ny, std::vector<double>(nx, 0.0));

    gamma11.assign(ny, std::vector<double>(nx, 0.0));
    gamma12.assign(ny, std::vector<double>(nx, 0.0));
    gamma22.assign(ny, std::vector<double>(nx, 0.0));

    norma11.assign(ny, std::vector<double>(nx, 0.0));
    norma12.assign(ny, std::vector<double>(nx, 0.0));

    tau11_aux.assign(ny, std::vector<double>(nx, 0.0));
    tau12_aux.assign(ny, std::vector<double>(nx, 0.0));
    tau22_aux.assign(ny, std::vector<double>(nx, 0.0));

    tau11_aux_prev.assign(ny, std::vector<double>(nx, 0.0));
    tau12_aux_prev.assign(ny, std::vector<double>(nx, 0.0));
    tau22_aux_prev.assign(ny, std::vector<double>(nx, 0.0));

    // Bingham allocations
    f11.assign(ny, std::vector<double>(nx, 0.0));
    f12.assign(ny, std::vector<double>(nx, 0.0));
    f22.assign(ny, std::vector<double>(nx, 0.0));

    tau11.assign(ny, std::vector<double>(nx, 0.0));
    tau12.assign(ny, std::vector<double>(nx, 0.0));
    tau22.assign(ny, std::vector<double>(nx, 0.0));

    gamma11.assign(ny, std::vector<double>(nx, 0.0));
    gamma12.assign(ny, std::vector<double>(nx, 0.0));
    gamma22.assign(ny, std::vector<double>(nx, 0.0));

    norma11.assign(ny, std::vector<double>(nx, 0.0));
    norma12.assign(ny, std::vector<double>(nx, 0.0));

    tau11_aux.assign(ny, std::vector<double>(nx, 0.0));
    tau12_aux.assign(ny, std::vector<double>(nx, 0.0));
    tau22_aux.assign(ny, std::vector<double>(nx, 0.0));

    tau11_aux_prev.assign(ny, std::vector<double>(nx, 0.0));
    tau12_aux_prev.assign(ny, std::vector<double>(nx, 0.0));
    tau22_aux_prev.assign(ny, std::vector<double>(nx, 0.0));

    // Derivatives and auxiliary tensors allocations
    du_dx.assign(ny, std::vector<double>(nx, 0.0));
    dv_dy.assign(ny, std::vector<double>(nx, 0.0));
    du_dy.assign(ny, std::vector<double>(nx, 0.0));
    dv_dx.assign(ny, std::vector<double>(nx, 0.0));
    D12.assign(ny, std::vector<double>(nx, 0.0));

    t11.assign(ny, std::vector<double>(nx, 0.0));
    t22.assign(ny, std::vector<double>(nx, 0.0));
    t12.assign(ny, std::vector<double>(nx, 0.0));

    // Build k-maps for u
    for (int j = 1; j <= ny_max; ++j)
        for (int i = 1; i <= nx_max; ++i)
            if (is_in_L_shape_u(j, i)) {
                u_node_to_k[{j, i}] = (int)u_k_to_node.size();
                u_k_to_node.push_back({j, i});
            }

    // Build k-maps for v
    for (int j = 1; j <= ny_max; ++j)
        for (int i = 1; i <= nx_max; ++i)
            if (is_in_L_shape_v(j, i)) {
                v_node_to_k[{j, i}] = (int)v_k_to_node.size();
                v_k_to_node.push_back({j, i});
            }
}

// ============================================================
// Precompute Matrices (done once)
// ============================================================

void LShapeSolver2D::precompute_matrices() {
    // ---- Matrix for U ----
    int Nu = (int)u_k_to_node.size();
    Au.resize(Nu, Nu); b_u.resize(Nu); sol_u.setZero(Nu);
    std::vector<Eigen::Triplet<double>> tripU;
    tripU.reserve(Nu * 5);

    for (int k = 0; k < Nu; ++k) {
        int j = u_k_to_node[k].first;
        int i = u_k_to_node[k].second;
        double diag = 4.0;

        // Left neighbour
        if (is_in_L_shape_u(j, i - 1)) {
            tripU.push_back(Eigen::Triplet<double>(k, u_node_to_k[{j, i - 1}], -1.0));
        } else {
            // Left wall i=1: Dirichlet u=0 -> diag=5
            diag += 1.0;
        }

        // Right neighbour
        if (is_in_L_shape_u(j, i + 1)) {
            tripU.push_back(Eigen::Triplet<double>(k, u_node_to_k[{j, i + 1}], -1.0));
        } else {
            if (j <= cfg.ny1 && i == cfg.nx2 - 1) {
                // Exit: Neumann u_x=0 -> diag=3
                diag -= 1.0;
            } else {
                // Wall at nx1 or step ceiling: u=0 -> diag=5
                diag += 1.0;
            }
        }

        // Bottom neighbour
        if (is_in_L_shape_u(j - 1, i)) {
            tripU.push_back(Eigen::Triplet<double>(k, u_node_to_k[{j - 1, i}], -1.0));
        } else {
            if (j == 1) {
                // Bottom: Symmetry u_y=0 -> diag=3
                diag -= 1.0;
            } else if (j == cfg.ny1 + 1 && i >= cfg.nx1) {
                // Step ceiling wide part: u=0 -> diag=5
                diag += 1.0;
            }
        }

        // Top neighbour
        if (is_in_L_shape_u(j + 1, i)) {
            tripU.push_back(Eigen::Triplet<double>(k, u_node_to_k[{j + 1, i}], -1.0));
        } else {
            // Top wall: u=0 -> diag=5
            diag += 1.0;
        }

        tripU.push_back(Eigen::Triplet<double>(k, k, diag));
    }
    Au.setFromTriplets(tripU.begin(), tripU.end());
    solver_u.compute(Au);

    // ---- Matrix for V ----
    int Nv = (int)v_k_to_node.size();
    Av.resize(Nv, Nv); b_v.resize(Nv); sol_v.setZero(Nv);
    std::vector<Eigen::Triplet<double>> tripV;
    tripV.reserve(Nv * 5);

    for (int k = 0; k < Nv; ++k) {
        int j = v_k_to_node[k].first;
        int i = v_k_to_node[k].second;
        double diag = 4.0;

        // Left neighbour
        if (is_in_L_shape_v(j, i - 1)) {
            tripV.push_back(Eigen::Triplet<double>(k, v_node_to_k[{j, i - 1}], -1.0));
        } else {
            if (i == 1) {
                // Left symmetry: v_x=0 -> diag=3
                diag -= 1.0;
            } else {
                // Left wall: v=0 -> diag=5
                diag += 1.0;
            }
        }

        // Right neighbour
        if (is_in_L_shape_v(j, i + 1)) {
            tripV.push_back(Eigen::Triplet<double>(k, v_node_to_k[{j, i + 1}], -1.0));
        } else {
            if (j < cfg.ny1 && i == cfg.nx2) {
                // Exit Neumann: v_x=0 -> diag=3
                diag -= 1.0;
            } else {
                // Wall at nx1: v=0 -> diag=5
                diag += 1.0;
            }
        }

        // Bottom neighbour
        if (is_in_L_shape_v(j - 1, i)) {
            tripV.push_back(Eigen::Triplet<double>(k, v_node_to_k[{j - 1, i}], -1.0));
        } else {
            // Bottom wall: v=0 -> diag=5
            diag += 1.0;
        }

        // Top neighbour
        if (is_in_L_shape_v(j + 1, i)) {
            tripV.push_back(Eigen::Triplet<double>(k, v_node_to_k[{j + 1, i}], -1.0));
        } else {
            // Top wall (inflow or no-slip): v=V_bc -> diag=5 (RHS gets 2*V_bc)
            diag += 1.0;
        }

        tripV.push_back(Eigen::Triplet<double>(k, k, diag));
    }
    Av.setFromTriplets(tripV.begin(), tripV.end());
    solver_v.compute(Av);
}

// ============================================================
// solve_poisson_u  (generic RHS)
// ============================================================

void LShapeSolver2D::solve_poisson_u(VectorField& u_field, const VectorField& rhs) {
    int Nu = (int)u_k_to_node.size();
    double h2 = cfg.h * cfg.h;
    for (int k = 0; k < Nu; ++k) {
        int j = u_k_to_node[k].first;
        int i = u_k_to_node[k].second;
        b_u(k) = rhs[j][i] * h2;
    }
    sol_u = solver_u.solve(b_u);
    for (int k = 0; k < Nu; ++k) {
        int j = u_k_to_node[k].first;
        int i = u_k_to_node[k].second;
        u_field[j][i] = sol_u(k);
    }
    // Ghost: right exit Neumann  u[j][nx2] = u[j][nx2-1]
    for (int j = 1; j <= cfg.ny1; ++j)
        u_field[j][cfg.nx2] = u_field[j][cfg.nx2 - 1];
}

// ============================================================
// solve_poisson_v  (with inflow BC at top)
// ============================================================

void LShapeSolver2D::solve_poisson_v(VectorField& v_field, const VectorField& rhs) {
    int Nv = (int)v_k_to_node.size();
    double h2 = cfg.h * cfg.h;
    for (int k = 0; k < Nv; ++k) {
        int j = v_k_to_node[k].first;
        int i = v_k_to_node[k].second;
        b_v(k) = rhs[j][i] * h2;
        // Dirichlet at top lid
        if (!is_in_L_shape_v(j + 1, i) && j == cfg.ny2 - 1) {
            double bc = (i <= cfg.lx_hole) ? cfg.v_inflow : 0.0;
            b_v(k) += 2.0 * bc;
        }
    }
    sol_v = solver_v.solve(b_v);
    for (int k = 0; k < Nv; ++k) {
        int j = v_k_to_node[k].first;
        int i = v_k_to_node[k].second;
        v_field[j][i] = sol_v(k);
    }
    // Ghost: top Dirichlet inflow  v[ny2][i] = v_inflow
    for (int i = 1; i <= cfg.nx1; ++i)
        v_field[cfg.ny2][i] = (i <= cfg.lx_hole) ? cfg.v_inflow : 0.0;
}

// ============================================================
// OPERATOR A
// ============================================================

void LShapeSolver2D::solve_poisson_u_A(VectorField& ut, const VectorField& p_field) {
    double h = cfg.h;
    int Nu = (int)u_k_to_node.size();
    for (int k = 0; k < Nu; ++k) {
        int j = u_k_to_node[k].first;
        int i = u_k_to_node[k].second;
        double pR = is_in_L_shape_p(j, i + 1) ? p_field[j][i + 1] : p_field[j][i];
        double pL = is_in_L_shape_p(j, i) ? p_field[j][i] : p_field[j][i + 1];
        b_u(k) = (pR - pL) * h;
    }
    sol_u = solver_u.solve(b_u);
    for (int k = 0; k < Nu; ++k) {
        int j = u_k_to_node[k].first;
        int i = u_k_to_node[k].second;
        ut[j][i] = sol_u(k);
    }
    for (int j = 1; j <= cfg.ny1; ++j)
        ut[j][cfg.nx2] = ut[j][cfg.nx2 - 1];
}

void LShapeSolver2D::solve_poisson_v_A(VectorField& vt, const VectorField& p_field) {
    double h = cfg.h;
    int Nv = (int)v_k_to_node.size();
    for (int k = 0; k < Nv; ++k) {
        int j = v_k_to_node[k].first;
        int i = v_k_to_node[k].second;
        double pT = is_in_L_shape_p(j + 1, i) ? p_field[j + 1][i] : p_field[j][i];
        double pB = is_in_L_shape_p(j, i) ? p_field[j][i] : p_field[j + 1][i];
        b_v(k) = (pT - pB) * h;
    }
    sol_v = solver_v.solve(b_v);
    for (int k = 0; k < Nv; ++k) {
        int j = v_k_to_node[k].first;
        int i = v_k_to_node[k].second;
        vt[j][i] = sol_v(k);
    }
}

void LShapeSolver2D::Operator_A(const VectorField& p_field) {
    zero_field(u_tilde);
    zero_field(v_tilde);
   

    // Dispatch u_tilde and v_tilde solves via pre-allocated pool
    auto fu_A = pool->enqueue([this, &p_field] { solve_poisson_u_A(u_tilde, p_field); });
    auto fv_A = pool->enqueue([this, &p_field] { solve_poisson_v_A(v_tilde, p_field); });
    fu_A.get();
    fv_A.get();

    // Ap = -div(u_tilde, v_tilde)  -- sequential (cheap loop)
    double op_h = 1.0/cfg.h;
    for (int j = 1; j <= ny_max; ++j) {
        for (int i = 1; i <= nx_max; ++i) {
            if (is_in_L_shape_p(j, i)) {
                // Apply boundary conditions explicitly for divergence
                double u_R = u_tilde[j][i];
                double u_L = u_tilde[j][i-1];
                double v_T = v_tilde[j][i];
                double v_B = v_tilde[j-1][i];

                // Right outflow Neumann BC  (u_nx2 = u_{nx2-1})
                if (j <= cfg.ny1 && i == cfg.nx2) u_R = u_tilde[j][i-1];
                
                // Top Dirichlet BC (homogeneous for operator A -> v = 0)
                if (j == cfg.ny2 && i <= cfg.nx1) v_T = 0.0;

                Ap[j][i] = -((u_R - u_L) * op_h + (v_T - v_B) * op_h);
            }
        }
    }
}

// ============================================================
// CG helpers
// ============================================================

void LShapeSolver2D::calculation_r() {
    // r1 = Ap - F  (Standard stationary residual for Uzawa)
    for (int j = 1; j <= ny_max; ++j)
        for (int i = 1; i <= nx_max; ++i)
            if (is_in_L_shape_p(j, i))
                r1[j][i] = Ap[j][i] - F[j][i];
}

// scalar_product_rr:
//   sum1 = (r1, r1),  sum2 = (r1, Ap)
//   sets r_cur = sum1, returns tau = sum1 / sum2   (like Caverna::scalar_product)
double LShapeSolver2D::scalar_product_rr() {
    double sum1 = 0.0, sum2 = 0.0, r = 0.0;
    for (int j = 1; j <= ny_max; ++j)
        for (int i = 1; i <= nx_max; ++i)
            if (is_in_L_shape_p(j, i)) {
                r = r1[j][i];
                sum1 += r * r;
                sum2 += r * Ap[j][i];
            }
    r_cur = sum1;
    return sum1 / sum2;  // = tau
}

// ============================================================
// solve_pressure_cg_custom
// ============================================================

void LShapeSolver2D::solve_pressure_cg_custom() {
    int max_iteration = 100;
    double tol = 1e-5;
    //std::cout << "CG Step 1:" << "\n";
    // Step 1: solve u_tilde, v_tilde via pre-allocated pool
    auto fut_u = pool->enqueue([this] { solve_poisson_u(u_tilde, fu); });
    auto fut_v = pool->enqueue([this] { solve_poisson_v(v_tilde, fv); });
    fut_u.get();
    fut_v.get();
    double prev , cur, r_prev_tmp, tau_prev_tmp;
   
    double op_h = 1.0/cfg.h;
    //std::cout << "CG Step 2:" << "\n";
    // Step 2: compute F[j][i] = -((u_R - u_L) / h + (v_T - v_B) / h) - g,  r1 = -F  (r_0)
    for (int j = 1; j <= ny_max; ++j) {
        for (int i = 1; i <= nx_max; ++i) {
            if (is_in_L_shape_p(j, i)) {
                // Apply boundary conditions explicitly for divergence
                double u_R = u_tilde[j][i];
                double u_L = u_tilde[j][i-1];
                double v_T = v_tilde[j][i];
                double v_B = v_tilde[j-1][i];

                double g = 0.0;

                // Right outflow Neumann BC  (u_nx2 = u_{nx2-1})
                if (j <= cfg.ny1 && i == cfg.nx2) u_R = u_tilde[j][i-1];
                
                // Top inflow Dirichlet BC
                if (j == cfg.ny2 && i <= cfg.nx1) {
                    v_T = (i <= cfg.lx_hole) ? cfg.v_inflow : 0.0;
                    // g = - V_boundary / h
                    //if (i <= cfg.lx_hole) g = -(cfg.v_inflow / h);
                }

                // F = -div - g
                F[j][i] = -((u_R - u_L) * op_h + (v_T - v_B) * op_h) ; //  - g
                r1[j][i] = -F[j][i];  // r_0 = -F
            }
        }
    }


    // Step 3: Ap = A*r0,  tau_0 = scalar_product_rr()
    Operator_A(r1);
    tau_cur = scalar_product_rr();  // also sets r_cur = (r0,r0)

    // p_0 = -tau_0 * r_0,  p_prev = 0
    for (int j = 1; j <= ny_max; ++j)
        for (int i = 1; i <= nx_max; ++i)
            if (is_in_L_shape_p(j, i)) {
                p[j][i]      = -tau_cur * r1[j][i];
                p_prev[j][i] = 0.0;
            }

    alpha_prev = 1.0;

    // Main CG loop (Chebyshev acceleration, like Caverna)
    for (int k = 1; k < max_iteration; ++k) {
        Operator_A(p);           // Ap = A p_k
        calculation_r();         // r1 = Ap - F

        double r_prev_tmp = r_cur;
        double tau_prev_tmp = tau_cur;

        Operator_A(r1);          // Ap = A r_k
        tau_cur = scalar_product_rr();  // sets r_cur = (r_k, r_k), returns new tau

        alpha_cur = 1.0 / (1.0 - (tau_cur * r_cur) / (tau_prev_tmp * r_prev_tmp * alpha_prev));

        if (k % 50 == 0) std::cout << "CG k=" << k << "  r=" << r_cur << "  tau=" << tau_cur << "  alpha=" << alpha_cur << "\n";

        for (int j = 1; j <= ny_max; ++j)
            for (int i = 1; i <= nx_max; ++i)
                if (is_in_L_shape_p(j, i)) {
                    double prev = p_prev[j][i]; double cur  = p[j][i];
                    p[j][i]= alpha_cur * cur + (1.0 - alpha_cur) * prev - alpha_cur * tau_cur * r1[j][i];
                    p_prev[j][i] = cur;
                }

        alpha_prev = alpha_cur;
        if (r_cur < tol) break;
    }
    if (r_cur > tol) std::cout << "CG did not converge: r=" << r_cur << "\n";

    // Step 4: Final velocities with pressure correction
    //   -Delta u = fu - Bx*p,   -Delta v = fv - By*p
    auto future_u = pool->enqueue([this] { solve_poisson_u_p(u, fu, p); });
    auto future_v = pool->enqueue([this] { solve_poisson_v_p(v, fv, p); });
    future_u.get();
    future_v.get();
}

// ============================================================
// solve_poisson_u_p:  -Delta u = fu - Bx*p
// ============================================================
void LShapeSolver2D::solve_poisson_u_p(VectorField& u_field, const VectorField& rhs, const VectorField& p_field) {
    double h = cfg.h;
    double h2 = h * h;
    int Nu = (int)u_k_to_node.size();
    for (int k = 0; k < Nu; ++k) {
        int j = u_k_to_node[k].first;
        int i = u_k_to_node[k].second;
        double pR = is_in_L_shape_p(j, i + 1) ? p_field[j][i + 1] : p_field[j][i];
        double pL = is_in_L_shape_p(j, i) ? p_field[j][i] : p_field[j][i + 1];
        b_u(k) = rhs[j][i] * h2 - (pR - pL) * h;
    }
    sol_u = solver_u.solve(b_u);
    for (int k = 0; k < Nu; ++k) {
        int j = u_k_to_node[k].first;
        int i = u_k_to_node[k].second;
        u_field[j][i] = sol_u(k);
    }
    for (int j = 1; j <= cfg.ny1; ++j)
        u_field[j][cfg.nx2] = u_field[j][cfg.nx2 - 1];
}

void LShapeSolver2D::solve_poisson_v_p(VectorField& v_field, const VectorField& rhs, const VectorField& p_field) {
    double h = cfg.h;
    double h2 = h * h;
    int Nv = (int)v_k_to_node.size();
    for (int k = 0; k < Nv; ++k) {
        int j = v_k_to_node[k].first;
        int i = v_k_to_node[k].second;
        double pT = is_in_L_shape_p(j + 1, i) ? p_field[j + 1][i] : p_field[j][i];
        double pB = is_in_L_shape_p(j, i) ? p_field[j][i] : p_field[j + 1][i];
        b_v(k) = rhs[j][i] * h2 - (pT - pB) * h;
        if (!is_in_L_shape_v(j + 1, i) && j == cfg.ny2 - 1) {
            double bc = (i <= cfg.lx_hole) ? cfg.v_inflow : 0.0;
            b_v(k) += 2.0 * bc;
        }
    }
    sol_v = solver_v.solve(b_v);
    for (int k = 0; k < Nv; ++k) {
        int j = v_k_to_node[k].first;
        int i = v_k_to_node[k].second;
        v_field[j][i] = sol_v(k);
    }
    for (int i = 1; i <= cfg.nx1; ++i)
        v_field[cfg.ny2][i] = (i <= cfg.lx_hole) ? cfg.v_inflow : 0.0;
}

// ============================================================
// ============================================================
// run_full_algorithm
// ============================================================

int LShapeSolver2D::run_full_algorithm(int max_iterations, double tolerance) {
    Timer timer;
    double t_sg;
    for (iter = 1; iter <= max_iterations; ++iter) {
        
        // Step 1: Pressure projection
        Timer timer_cg_custom;
        solve_pressure_cg_custom();
        t_sg = timer_cg_custom.elapsed();
        // Step 2: Bingham updates
        update_gamma();
        double current_diff = update_tau();
        diff_max.push_back(current_diff);
        update_f();

        std::cout << "Iteration " << iter << " | Diff Max: " << current_diff << " | Time CG:" << t_sg << " s " << " | Time Algorithm: " << timer.elapsed() << "s" << std::endl;

        // Periodic save for plotting
        if (iter % 1 == 0) {
            //norma_11_12(tau11, tau12, tau22);
            std::string prefix = "iter_" + std::to_string(iter) + "_";
            save_field(prefix + "norma11.txt", norma11, 5);
            save_field(prefix + "norma12.txt", norma12, 4);
            /*save_field(prefix + "p.txt", p, 3);
            save_field(prefix + "u.txt", u, 1);
            save_field(prefix + "v.txt", v, 2);*/
           
             std::ofstream diff_file(out_path("diff_max.txt"));
            if (diff_file.is_open()) {
                for (double d : diff_max) diff_file << d << "\n";
                diff_file.close();
            }
        }    
        
        if (current_diff < tolerance) { 
            std::cout << "Converged in " << iter << " iterations." << std::endl; 
            break; 
        }
    }

    std::string prefix = "iter_" + std::to_string(iter) + "_";
    save_field(prefix + "norma11.txt", norma11, 5);
    save_field(prefix + "norma12.txt", norma12, 4);
    
    std::ofstream diff_file(out_path("diff_max.txt"));
    if (diff_file.is_open()) {
        for (double d : diff_max) diff_file << d << "\n";
        diff_file.close();
    }
    return 0;
}

// ============================================================
// Bingham Steps
// ============================================================

double LShapeSolver2D::avg4(const VectorField& field, int j, int i) {
    auto safe_val = [&](int cj, int ci) {
        if (cj < 0) cj = 0;
        if (cj >= field.size()) cj = field.size() - 1;
        if (ci < 0) ci = 0;
        if (ci >= field[0].size()) ci = field[0].size() - 1;
        return field[cj][ci];
    };
    return 0.25 * (safe_val(j, i) + safe_val(j, i-1) + safe_val(j-1, i) + safe_val(j-1, i-1));
}

double LShapeSolver2D::avg2x(const VectorField& field, int j, int i) {
    auto safe_val = [&](int cj, int ci) {
        if (cj < 0) cj = 0;
        if (cj >= field.size()) cj = field.size() - 1;
        if (ci < 0) ci = 0;
        if (ci >= field[0].size()) ci = field[0].size() - 1;
        return field[cj][ci];
    };
    return 0.5 * (safe_val(j, i) + safe_val(j, i-1));
}

double LShapeSolver2D::avg2y(const VectorField& field, int j, int i) {
    auto safe_val = [&](int cj, int ci) {
        if (cj < 0) cj = 0;
        if (cj >= field.size()) cj = field.size() - 1;
        if (ci < 0) ci = 0;
        if (ci >= field[0].size()) ci = field[0].size() - 1;
        return field[cj][ci];
    };
    return 0.5 * (safe_val(j, i) + safe_val(j-1, i));
}

void LShapeSolver2D::norma_11_12(const VectorField& q11, const VectorField& q12, const VectorField& q22)
{
    double T11, T22, T12;
    int nx = nx_max;
    int ny = ny_max;

    // norma 11
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            T11 = q11[j][i];  T22 = q22[j][i];  T12 = avg4(q12, j, i);
            norma11[j][i] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);
        }
    }

    for (int j = 1; j <= ny; ++j) {
        T11 = q11[j][0];   T22 = q22[j][0]; T12 = avg2y(q12, j, 0); 
        norma11[j][0] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);

        T11 = q11[j][nx + 1];   T22 = q22[j][nx + 1]; T12 = avg2y(q12, j, nx);
        norma11[j][nx + 1] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);
    }

    for (int i = 1; i <= nx; ++i) {
        T11 = q11[0][i];   T22 = q22[0][i]; T12 = avg2x(q12, 0, i);
        norma11[0][i] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);
        
        T11 = q11[ny + 1][i];   T22 = q22[ny + 1][i]; T12 = avg2x(q12, ny, i);
        norma11[ny + 1][i] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);
    }
    
    // Corners for norma 11
    norma11[0][0] = 0.0; norma11[0][nx+1] = 0.0;
    norma11[ny+1][0] = 0.0; norma11[ny+1][nx+1] = 0.0;

    // norma 12
    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            T11 = avg4(q11, j+1, i+1);  T22 = avg4(q22, j+1, i+1); T12 = q12[j][i];
            norma12[j][i] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);
        }
    }

    for (int j = 1; j < ny; ++j) {
        T11 = avg2y(q11, j+1, 0);  T22 = avg2y(q22, j+1, 0); T12 = q12[j][0];
        norma12[j][0] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);

        T11 = avg2y(q11, j+1, nx+1);  T22 = avg2y(q22, j+1, nx+1); T12 = q12[j][nx];
        norma12[j][nx] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);
    }

    for (int i = 1; i < nx; ++i) {
        T11 = avg2x(q11, 0, i+1);  T22 = avg2x(q22, 0, i+1); T12 = q12[0][i];
        norma12[0][i] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);

        T11 = avg2x(q11, ny+1, i+1);  T22 = avg2x(q22, ny+1, i+1); T12 = q12[ny][i];
        norma12[ny][i] = std::sqrt(T11*T11 + T22*T22 + 2.0*T12*T12);
    }
}

void LShapeSolver2D::update_gamma() {
    // Step 1: Compute derivatives
    compute_derivatives();
    double dom = 2.0 * (r_param + mu);
    
    // Step 2: Update auxiliary tensors t_ij
    update_t();

    // Step 3: Compute norms for gamma updates
    norma_11_12(t11, t12, t22);

    // Step 4: Update gamma12 (nodes)
    for (int j = 0; j <= ny_max; ++j) {
        for (int i = 0; i <= nx_max; ++i) {
            update_gamma12(norma12[j][i], dom, i, j);
        }
    }

    // Step 5: Update gamma11, gamma22 (centers)
    for (int j = 0; j <= ny_max + 1; ++j) {
        for (int i = 0; i <= nx_max + 1; ++i) {
            update_gamma11_22(norma11[j][i], dom, i, j);
        }
    }
}

void LShapeSolver2D::compute_derivatives() {
    double h = cfg.h;
    // 1. du_dx and dv_dy (Centers)
    for (int j = 1; j <= ny_max; ++j) {
        for (int i = 1; i <= nx_max; ++i) {
            if (is_in_L_shape_p(j, i)) {
                du_dx[j][i] = (u[j][i] - u[j][i - 1]) / h;
                dv_dy[j][i] = (v[j][i] - v[j - 1][i]) / h;
            } else {
                du_dx[j][i] = 0.0; dv_dy[j][i] = 0.0;
            }
        }
    }

    // Boundary stencils for du_dx, dv_dy
    // Left wall i=0 (No-slip u=0, v=0)
    for (int j = 1; j <= ny_max; ++j) {
        if (is_in_L_shape_p(j, 1)) {
            du_dx[j][0] = (-3.0 * 0.0 + 4.0 * u[j][1] - u[j][2]) / (2.0 * h);
        }
    }
    // Bottom wall j=0 (Symmetry u_y=0, v=0)
    for (int i = 1; i <= nx_max; ++i) {
        if (is_in_L_shape_p(1, i)) {
            dv_dy[0][i] = (-3.0 * 0.0 + 4.0 * v[1][i] - v[2][i]) / (2.0 * h);
        }
    }
    // Top wall j=ny2 (Inflow u=0, v=v_inflow)
    for (int i = 1; i <= cfg.nx1; ++i) {
        double vy_top = (i <= cfg.lx_hole) ? cfg.v_inflow : 0.0;
        dv_dy[ny_max + 1][i] = (3.0 * vy_top - 4.0 * v[ny_max - 1][i] + v[ny_max - 2][i]) / (2.0 * h);
    }
    
    // Step internal boundaries (i=nx1 for j>ny1, j=ny1 for i>nx1)
    for (int j = cfg.ny1 + 1; j <= ny_max; ++j) {
        // Vertical step wall: du/dx at x=nx1
        du_dx[j][cfg.nx1] = (3.0 * 0.0 - 4.0 * u[j][cfg.nx1 - 1] + u[j][cfg.nx1 - 2]) / (2.0 * h);
    }
    for (int i = cfg.nx1 + 1; i <= cfg.nx2; ++i) {
        // Horizontal step ceiling: dv/dy at y=ny1
        dv_dy[cfg.ny1][i] = (3.0 * 0.0 - 4.0 * v[cfg.ny1 - 1][i] + v[cfg.ny1 - 2][i]) / (2.0 * h);
    }

    // Right exit i=nx2 (For now set to 0 as requested)
    for (int j = 1; j <= cfg.ny1; ++j) {
        du_dx[j][cfg.nx2] = 0.0;
        dv_dy[j][cfg.nx2] = 0.0;
    }

    // 2. du_dy and dv_dx (Nodes)
    for (int j = 1; j <= ny_max - 1; ++j) {
        for (int i = 1; i <= nx_max - 1; ++i) {
            bool has_u = is_in_L_shape_u(j, i) && is_in_L_shape_u(j + 1, i);
            bool has_v = is_in_L_shape_v(j, i) && is_in_L_shape_v(j, i + 1);
            if (has_u && has_v) {
                du_dy[j][i] = (u[j + 1][i] - u[j][i]) / h;
                dv_dx[j][i] = (v[j][i + 1] - v[j][i]) / h;
            } else {
                du_dy[j][i] = 0.0; dv_dx[j][i] = 0.0;
            }
        }
    }

    // Wall stencils for nodes (8/9/3 formula)
    // Horizontal walls (Symmetry bottom j=0, No-slip top j=ny2, Step ceiling j=ny1)
    for (int i = 1; i < nx_max; ++i) {
        // Bottom symmetry: du/dy = 0
        if (is_in_L_shape_u(1, i)) du_dy[0][i] = 0.0;
        
        // Top lid (no-slip u=0)
        if (i <= cfg.nx1 && is_in_L_shape_u(ny_max, i)) {
            du_dy[ny_max][i] = (8.0 * 0.0 - 9.0 * u[ny_max][i] + (ny_max > 1 ? u[ny_max - 1][i] : 0.0)) / (3.0 * h);
        }
        
        // Step ceiling (no-slip u=0)
        if (i > cfg.nx1 && i <= cfg.nx2 && is_in_L_shape_u(cfg.ny1, i)) {
             du_dy[cfg.ny1][i] = (8.0 * 0.0 - 9.0 * u[cfg.ny1][i] + u[cfg.ny1-1][i]) / (3.0 * h);
        }
    }

    // Vertical walls (No-slip left i=0, Step side i=nx1)
    for (int j = 1; j < ny_max; ++j) {
        // Left wall (v=0)
        if (is_in_L_shape_v(j, 1)) {
            dv_dx[j][0] = (-8.0 * 0.0 + 9.0 * v[j][1] - v[j][2]) / (3.0 * h);
        }
        // Step side (v=0)
        if (j > cfg.ny1 && is_in_L_shape_v(j, cfg.nx1)) {
            dv_dx[j][cfg.nx1] = (8.0 * 0.0 - 9.0 * v[j][cfg.nx1] + v[j][cfg.nx1-1]) / (3.0 * h);
        }
        // Right exit (i=nx2) - set to 0 for now as requested
        if (j <= cfg.ny1) {
            dv_dx[j][cfg.nx2] = 0.0;
        }
    }
}

void LShapeSolver2D::update_t() {
    // Centers
    for (int j = 0; j <= ny_max + 1; ++j) {
        for (int i = 0; i <= nx_max + 1; ++i) {
            t11[j][i] = tau11[j][i] + 2.0 * r_param * du_dx[j][i];
            t22[j][i] = tau22[j][i] + 2.0 * r_param * dv_dy[j][i];
        }
    }
    // Nodes
    for (int j = 0; j <= ny_max; ++j) {
        for (int i = 0; i <= nx_max; ++i) {
            t12[j][i] = tau12[j][i] + r_param * (du_dy[j][i] + dv_dx[j][i]);
            D12[j][i] = (du_dy[j][i] + dv_dx[j][i]) / 2.0;
        }
    }
}

void LShapeSolver2D::update_gamma12(double modul, double dom, int i, int j) {
    if (modul < tau_s)
        gamma12[j][i] = 0;
    else
        gamma12[j][i] = (1.0 - tau_s / modul) * t12[j][i] / dom;
}

void LShapeSolver2D::update_gamma11_22(double modul, double dom, int i, int j) {
    if (modul < tau_s) {
        gamma11[j][i] = 0; gamma22[j][i] = 0;
    } else {
        gamma11[j][i] = (1.0 - tau_s / modul) * t11[j][i] / dom;
        gamma22[j][i] = (1.0 - tau_s / modul) * t22[j][i] / dom;
    }
}

double LShapeSolver2D::update_tau() {
    double diff = 0.0;
    double T11, T22, T12;

    // Step 1: Update tau11, tau22 and compute diff (Centers)
    for (int j = 0; j <= ny_max + 1; ++j) {
        for (int i = 0; i <= nx_max + 1; ++i) {
            T11 = (du_dx[j][i] - gamma11[j][i]);
            T22 = (dv_dy[j][i] - gamma22[j][i]);

            tau11[j][i] += 2.0 * rho * T11;
            tau22[j][i] += 2.0 * rho * T22;

            diff = std::max({ diff, std::abs(T11), std::abs(T22) });
        }
    }

    // Step 2: Update tau12 and compute diff (Nodes)
    for (int j = 0; j <= ny_max; ++j) {
        for (int i = 0; i <= nx_max; ++i) {
            T12 = (D12[j][i] - gamma12[j][i]);
            tau12[j][i] += 2.0 * rho * T12;

            diff = std::max(diff, std::abs(T12));
        }
    }

    // Step 3: Compute norms to see true tau state
    norma_11_12(tau11, tau12, tau22);
    
    return diff;
}
void LShapeSolver2D::zero_field(VectorField& field) {
    for (auto& row : field) {
        std::fill(row.begin(), row.end(), 0.0);
    }
}

void LShapeSolver2D::update_f() {
    double inv_h = 1.0 / cfg.h;
    double inv_dom = 1.0 / (r_param); // Масштабный множитель уравнения

    // Step 1: Calculate Q = tau - 2 * r * gamma and store in temporary buffer t
    // Нам нужны значения в фиктивных слоях для корректных разностей на границах
    for (int j = 0; j <= ny_max + 1; ++j) {
        for (int i = 0; i <= nx_max + 1; ++i) {
            t11[j][i] = tau11[j][i] - 2.0 * r_param * gamma11[j][i];
            t22[j][i] = tau22[j][i] - 2.0 * r_param * gamma22[j][i];
        }
    }
    for (int j = 0; j <= ny_max; ++j) {
        for (int i = 0; i <= nx_max; ++i) {
            t12[j][i] = tau12[j][i] - 2.0 * r_param * gamma12[j][i];
        }
    }

    // Экстраполяция t12 на внешние границы (i=0, i=nx, j=0, j=ny) для устойчивости
    for (int j = 0; j <= ny_max; ++j) {
        t12[j][0] = t12[j][1];
        if (cfg.nx2 <= nx_max) t12[j][cfg.nx2] = t12[j][cfg.nx2 - 1];
    }
    for (int i = 0; i <= nx_max; ++i) {
        t12[0][i] = t12[1][i];
        if (cfg.ny2 <= ny_max) t12[cfg.ny2][i] = t12[cfg.ny2 - 1][i];
    }

    // Step 2: Compute fu and fv as divergence of Q
    for (int j = 1; j <= ny_max; ++j) {
        for (int i = 1; i <= nx_max; ++i) {

            // fu (в узлах скорости u)
            if (is_in_L_shape_u(j, i)) {
                // Дивергенция в узле u_ji: dQ11/dx + dQ12/dy
                double dQ11_dx = (t11[j][i + 1] - t11[j][i]) * inv_h;
                double dQ12_dy = (t12[j][i] - t12[j - 1][i]) * inv_h;
                fu[j][i] = (dQ11_dx + dQ12_dy);
            }
            else {
                fu[j][i] = 0.0;
            }

            // fv (в узлах скорости v)
            if (is_in_L_shape_v(j, i)) {
                // Дивергенция в узле v_ji: dQ22/dy + dQ12/dx
                double dQ22_dy = (t22[j + 1][i] - t22[j][i]) * inv_h;
                double dQ12_dx = (t12[j][i] - t12[j][i - 1]) * inv_h;
                fv[j][i] = (dQ22_dy + dQ12_dx);
            }
            else {
                fv[j][i] = 0.0;
            }
        }
    }
}


// ============================================================
// save_field
// ============================================================

void LShapeSolver2D::save_field(const std::string& filename, const VectorField& data1, int boll) {
    CreateDirectoryA(cfg.output_dir.c_str(), NULL);
    std::string path = out_path(filename);
    std::ofstream file(path);
    if (!file.is_open()) return;
    file << std::fixed << std::setprecision(8);
    double h  = cfg.h;
    double h2 = h / 2.0;

    if (boll == 1) { // u
        for (int j = 1; j <= cfg.ny2; ++j) {
            int cur_nx = (j <= cfg.ny1) ? cfg.nx2 : cfg.nx1;
            for (int i = 0; i <= cur_nx; ++i)
                file << (i * h) << " " << (j * h - h2) << " " << data1[j][i] << "\n";
            file << "\n";
        }
    } else if (boll == 2) { // v
        for (int j = 0; j <= cfg.ny2; ++j) {
            int cur_nx = (j < cfg.ny1) ? cfg.nx2 : cfg.nx1;
            for (int i = 1; i <= cur_nx; ++i)
                file << (i * h - h2) << " " << (j * h) << " " << data1[j][i] << "\n";
            file << "\n";
        }
    } else if (boll == 3) { // p
        for (int j = 1; j <= cfg.ny2; ++j) {
            int cur_nx = (j <= cfg.ny1) ? cfg.nx2 : cfg.nx1;
            for (int i = 1; i <= cur_nx; ++i)
                file << (i * h - h2) << " " << (j * h - h2) << " " << data1[j][i] << "\n";
            file << "\n";
        }
    } else if (boll == 4) { // nodes: x=i*h, y=j*h
        for (int j = 0; j <= cfg.ny2; ++j) {
            int cur_nx = (j <= cfg.ny1) ? cfg.nx2 : cfg.nx1;
            for (int i = 0; i <= cur_nx; ++i)
                file << (i * h) << " " << (j * h) << " " << data1[j][i] << "\n";
            file << "\n";
        }
    } else if (boll == 5) {
        int nx = nx_max;
        int ny = ny_max;
        double Lx1 = nx * h;
        double Ly1 = ny * h;
        double hx = h;
        double hy = h;
        
        file << 0.0 << " " << 0.0 << " " << data1[0][0] << std::endl;
        for (int i = 1; i <= nx; ++i) {
            double x = i * hx - hx / 2.0; double y = 0;
            file << x << " " << y << " " << data1[0][i] << std::endl;
        }
        file << Lx1 << " " << 0.0 << " " << data1[0][nx + 1] << std::endl;
        for (int j = 1; j <= ny; ++j) {
            double x = 0; double y = j * hy - hy / 2.0;
            file << x << " " << y << " " << data1[j][0] << std::endl;

            for (int i = 1; i <= nx; ++i) {
                x = i * hx - hx / 2.0; y = j * hy - hy / 2.0;
                // Leave zeros directly for missing L_shape parts just in case
                file << x << " " << y << " " << data1[j][i] << std::endl;
            }
            double xb = Lx1; y = j * hy - hy / 2.0;
            file << xb << " " << y << " " << data1[j][nx + 1] << std::endl;
            file << std::endl;
        }
        file << 0.0 << " " << Ly1 << " " << data1[ny + 1][0] << std::endl;
        for (int i = 1; i <= nx; ++i) {
            double x = i * hx - hx / 2.0; double y = Ly1;
            file << x << " " << y << " " << data1[ny+1][i] << std::endl;
        }
        // Corners
        file << Lx1 << " " << Ly1 << " " << data1[ny + 1][nx + 1] << std::endl;
    }
    file.close();
}
