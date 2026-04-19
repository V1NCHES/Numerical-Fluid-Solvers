#include "CavernaSolver2D.h"

// ///////////////////////////////////////////////////////////////////////////

void CavernaSolver2D::solve_pressure_cg_custom()
{
    double cur, prev;
    int max_iteration = 100;
    double tol = 1e-12;
    int k = 1;

    auto future_u_tilde = pool->enqueue([this] { solve_poisson_u_caverna(u_tilde, fu); });
    auto future_v_tilde = pool->enqueue([this] { solve_poisson_v(v_tilde, fv); });
    future_u_tilde.get();
    future_v_tilde.get();
    save_all(1);

    // Step 2: F
    for (int j = 1; j <= ny; ++j)
        for (int i = 1; i <= nx; ++i)
        {
            F[j][i] = -( (u_tilde[j][i] - u_tilde[j][i - 1]) / hx + (v_tilde[j][i] - v_tilde[j - 1][i]) / hy); r1[j][i] = -F[j][i]; // -g[j][i]
        } // r_0   
    Operator_A(r1); // Ar_0 
    tau_cur = scalar_product();
    //std::cout << "\nk = 0; r_cur = " << r_cur << " ; tau_cur = " << tau_cur << " sqrt(r_cur) = " << sqrt(r_cur) << "\n";
    for (int j = 1; j < ny + 1; ++j)
        for (int i = 1; i < nx + 1; ++i)
            p[j][i] = -tau_cur * r1[j][i], p_prev[j][i] = 0.0; // p_0 - tau_1 * r_0               

    alpha_prev = 1.0;
    // Step 3: Ищем p_k: A p_k = F
    for (k = 1; k < max_iteration; ++k)
    {
        Operator_A(p); //A p_k
        calculation_r(); // r_k
        r_prev = r_cur; tau_prev = tau_cur; // r_k-1 = r_k tau_k-1 = tau_k
        Operator_A(r1); //Ar_k
        tau_cur = scalar_product(); //r_k, tau_k
        alpha_cur = 1.0 / (1.0 - (tau_cur * r_cur) / (tau_prev * r_prev * alpha_prev)); // alpha_k+1
        if (k % 50 == 0) std::cout << "k =  " << k << "; r_cur = " << r_cur << " ; tau_cur = " << tau_cur << " ; alpha_cur = " << alpha_cur << " sqrt(r_cur) = " << sqrt(r_cur) << " ; k = " << k << "\n";
        for (int j = 1; j < ny + 1; ++j)
            for (int i = 1; i < nx + 1; ++i)
            {
                prev = p_prev[j][i]; cur = p[j][i];
                p[j][i] = alpha_cur * cur + (1 - alpha_cur) * prev - alpha_cur * tau_cur * r1[j][i];
                p_prev[j][i] = cur;
            }
        alpha_prev = alpha_cur; // alpha_k = alpha_k+1
        if (r_cur < tol) { break; } // std::cout << "final: k = " << k << " ; r_cur = " << r_cur << "\n";
    }
    if (r_cur > tol) { std::cout << "Error: r_cur = " << r_cur << "\n"; }
    //if (sqrt(r_cur) > 0.0000001) { std::cout << "k =  " << k << "; r_cur = " << r_cur << " ; tau_cur = " << tau_cur << " ; alpha_cur = " << alpha_cur << " sqrt(r_cur) = " << sqrt(r_cur) << " ; k = " << k << "\n"; }


    //antysymetric_x(p, nx + 1, ny + 1);
    // Step 4: Final velocities (incorporating pressure gradient)     
    auto future_u = pool->enqueue([this] { solve_poisson_u_p(u, fu, p); });
    auto future_v = pool->enqueue([this] { solve_poisson_v_p(v, fv, p); });
    future_u.get();
    future_v.get();

    // Step 5: Stream function
    /*for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            f_psi[j][i] = (v[j][i + 1] - v[j][i]) / hx - (u[j + 1][i] - u[j][i]) / hy;
        }
    }
    solve_poisson_psi(psi, f_psi);*/
   /* if (iter == 1)
    {
        save_all(1);
    }*/
}


double CavernaSolver2D::scalar_product()
{
    double sum1 = 0.0, sum2 = 0, r = 0;
    for (int j = 1; j < ny + 1; ++j)
        for (int i = 1; i < nx + 1; ++i)
        {
            r = r1[j][i];
            sum1 += r * r;
            sum2 += r * Ap[j][i];
        }
    r_cur = sum1;
    //std::cout << "\nscalar_product:  sum1 = " << sum1 << " ; sum2 = " << sum2 << "\n";
    return sum1 / sum2;
}

void CavernaSolver2D::calculation_r()
{
    for (int j = 1; j < ny + 1; ++j)
        for (int i = 1; i < nx + 1; ++i)
            r1[j][i] = Ap[j][i] - F[j][i];
}

/*
for (int j = 1; j <= ny; ++j) {
            for (int i = 1; i <= nx - 1; ++i) {
                Bx_p[j][i] = (p_field[j][i + 1] - p_field[j][i]) / hx;
            }
        }

        for (int j = 1; j <= ny - 1; ++j) {
            for (int i = 1; i <= nx; ++i) {
                By_p[j][i] = (p_field[j + 1][i] - p_field[j][i]) / hy;
            }
        }
*/

void CavernaSolver2D::Operator_A(const VectorField& p_field) {
    Timer timer_A;
    // Запускаем расчет Bx и решение по U в первом потоке
    auto future_u = pool->enqueue([this, &p_field] {solve_poisson_u_A(u_tilde, p_field);});
    // Запускаем расчет By и решение по V во втором потоке
    auto future_v = pool->enqueue([this, &p_field] {solve_poisson_v_A(v_tilde, p_field);});
    // Ожидаем завершения обеих задач
    future_u.get();
    future_v.get();

    for (int j = 1; j <= ny; ++j)
    {
        for (int i = 1; i <= nx; ++i) {
            //Ap[j][i] = -(BxT[j][i] + ByT[j][i]);
            Ap[j][i] = -((u_tilde[j][i] - u_tilde[j][i - 1]) / hx + (v_tilde[j][i] - v_tilde[j - 1][i]) / hy);
        }
    }

    //std::cout << " | Time  A(): " << timer_A.elapsed() << "s" << std::endl;
}




// ///////////////////////////////////////////////////////////////////////////


void CavernaSolver2D::calculation_g()
{
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            g[j][i] = 0;//-(u_tilde[j][i] - u_tilde[j][i - 1]) / hx - (v_tilde[j][i] - v_tilde[j - 1][i]) / hy;
        }
    }
}

void CavernaSolver2D::precompute_matrices() {
    // 1. Matrix for U: (nx-1) x ny
    int Nu = (nx - 1) * (ny);

    sol_u.setZero(Nu); // Теперь буфер памяти внутри sol_u готов и имеет нужный размер

    Au.resize(Nu, Nu); b_u.resize(Nu);
    std::vector<Eigen::Triplet<double>> tripletsU;
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            double diag = 4.0;
            if (j == 1) diag += 1.0;  // u_{i,0} = -u_{i,1} -> -3 total in Y part? No, -1 extra.
            if (j == ny) diag += 1.0; // u_{i,ny+1} = 2u_top - u_{i,ny}
            tripletsU.push_back(Eigen::Triplet<double>(k, k, diag));
            if (i > 1) tripletsU.push_back(Eigen::Triplet<double>(k, k - 1, -1.0));
            if (i < nx - 1) tripletsU.push_back(Eigen::Triplet<double>(k, k + 1, -1.0));
            if (j > 1) tripletsU.push_back(Eigen::Triplet<double>(k, k - (nx - 1), -1.0));
            if (j < ny) tripletsU.push_back(Eigen::Triplet<double>(k, k + (nx - 1), -1.0));
        }
    }
    Au.setFromTriplets(tripletsU.begin(), tripletsU.end());
    solver_u.compute(Au);

    // 2. Matrix for V: nx x (ny-1)
    int Nv = (nx) * (ny - 1);
    sol_v.setZero(Nv);

    Av.resize(Nv, Nv); b_v.resize(Nv);
    std::vector<Eigen::Triplet<double>> tripletsV;
    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx; ++i) {
            int k = (j - 1) * nx + (i - 1);
            double diag = 4.0;
            if (i == 1) diag += 1.0;  // v_{0,j} = -v_{1,j}
            if (i == nx) diag += 1.0; // v_{nx+1,j} = -v_{nx,j}
            tripletsV.push_back(Eigen::Triplet<double>(k, k, diag));
            if (i > 1) tripletsV.push_back(Eigen::Triplet<double>(k, k - 1, -1.0));
            if (i < nx) tripletsV.push_back(Eigen::Triplet<double>(k, k + 1, -1.0));
            if (j > 1) tripletsV.push_back(Eigen::Triplet<double>(k, k - nx, -1.0));
            if (j < ny - 1) tripletsV.push_back(Eigen::Triplet<double>(k, k + nx, -1.0));
        }
    }
    Av.setFromTriplets(tripletsV.begin(), tripletsV.end());
    solver_v.compute(Av);

    // 3. Matrix for PSI: (nx-1) x (ny-1)
    int Npsi = (nx - 1) * (ny - 1);
    sol_psi.setZero(Npsi);
    Apsi.resize(Npsi, Npsi); b_psi.resize(Npsi);
    std::vector<Eigen::Triplet<double>> tripletsPsi;
    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            tripletsPsi.push_back(Eigen::Triplet<double>(k, k, 4.0));
            if (i > 1) tripletsPsi.push_back(Eigen::Triplet<double>(k, k - 1, -1.0));
            if (i < nx - 1) tripletsPsi.push_back(Eigen::Triplet<double>(k, k + 1, -1.0));
            if (j > 1) tripletsPsi.push_back(Eigen::Triplet<double>(k, k - (nx - 1), -1.0));
            if (j < ny - 1) tripletsPsi.push_back(Eigen::Triplet<double>(k, k + (nx - 1), -1.0));
        }
    }
    Apsi.setFromTriplets(tripletsPsi.begin(), tripletsPsi.end());
    solver_psi.compute(Apsi);
}
// ///////////////////////////////////////////////////////////////////////////
void CavernaSolver2D::solve_poisson_u(VectorField& u_field, const VectorField& rhs_field) {
    int Nu = (nx - 1) * (ny);

    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            b_u(k) = rhs_field[j][i] * h2; //if (j == ny) {b_u(k) -= 2.0; // Boundary u=1 at the top: Adds 2/h^2 to Laplacian side -> subtract 2 from RHS vector k//}
        }
    }

    sol_u = solver_u.solve(b_u);
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            u_field[j][i] = sol_u(k);
        }
    }
}

void CavernaSolver2D::solve_poisson_u_caverna(VectorField& u_field, const VectorField& rhs_field) {
    int Nu = (nx - 1) * (ny);

    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            b_u(k) = rhs_field[j][i] * h2; //if (j == ny) {//    b_u(k) += 2.0; // Boundary u=1 at the top: Adds 2/h^2 to Laplacian side -> subtract 2 from RHS vector k //}
        }
    }
    int j = ny;
    for (int i = 1; i <= nx - 1; ++i)
    {
        int k = (j - 1) * (nx - 1) + (i - 1); //b_u(k) = 2.0;
        b_u(k) = rhs_field[j][i] * h2 + 2.0; // Boundary u=1 at the top: Adds 2/h^2 to Laplacian side -> subtract 2 from RHS vector k
    }

    sol_u = solver_u.solve(b_u);
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            u_field[j][i] = sol_u(k);
        }
    }
}

void CavernaSolver2D::solve_poisson_v(VectorField& v_field, const VectorField& rhs_field) {
    int Nv = (nx) * (ny - 1);

    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx; ++i) {
            int k = (j - 1) * nx + (i - 1);
            b_v(k) = rhs_field[j][i] * h2; // Boundaries v=0 on all walls: no adjustment to b needed (implicit 0)
        }
    }

    sol_v = solver_v.solve(b_v);
    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx; ++i) {
            int k = (j - 1) * nx + (i - 1);
            v_field[j][i] = sol_v(k);
        }
    }
}

void CavernaSolver2D::solve_poisson_psi(VectorField& psi_field, const VectorField& rhs_field) {
    int Npsi = (nx - 1) * (ny - 1);

    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            b_psi(k) = rhs_field[j][i] * h2;
        }
    }
    sol_psi = solver_psi.solve(b_psi);
    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            psi_field[j][i] = sol_psi(k);
        }
    }
}



void CavernaSolver2D::solve_poisson_u_p(VectorField& u_field, const VectorField& rhs_field, VectorField& p_field) {
    int Nu = (nx - 1) * ny;

    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            b_u(k) = (-(p_field[j][i + 1] - p_field[j][i]) * hx + rhs_field[j][i] * h2); //if (j == ny) {//    b_u(k) += 2.0; // Boundary u=1 at the top: Adds 2/h^2 to Laplacian side -> subtract 2 from RHS vector k //}
        }
    }
    int j = ny;
    for (int i = 1; i <= nx - 1; ++i)
    {
        int k = (j - 1) * (nx - 1) + (i - 1); //b_u(k) = 2.0;
        b_u(k) = (-(p_field[j][i + 1] - p_field[j][i]) * hx + rhs_field[j][i] * h2) + 2.0;  // Boundary u=1 at the top: Adds 2/h^2 to Laplacian side -> subtract 2 from RHS vector k
    }

    sol_u = solver_u.solve(b_u);
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            u_field[j][i] = sol_u(k);
        }
    }
}

void CavernaSolver2D::solve_poisson_v_p(VectorField& v_field, const VectorField& rhs_field, VectorField& p_field) {
    int Nv = nx * (ny - 1);

    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx; ++i) {
            int k = (j - 1) * nx + (i - 1);
            b_v(k) = (-(p_field[j + 1][i] - p_field[j][i]) * hx + rhs_field[j][i] * h2);// Boundaries v=0 on all walls: no adjustment to b needed (implicit 0)
        }
    }
    sol_v = solver_v.solve(b_v);
    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx; ++i) {
            int k = (j - 1) * nx + (i - 1);
            v_field[j][i] = sol_v(k);
        }
    }
}

void CavernaSolver2D::solve_poisson_u_A(VectorField& u_field, const VectorField& P) {
    int Nu = (nx - 1) * (ny);

    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            b_u(k) = (P[j][i + 1] - P[j][i]) * hx;
        }
    }
    sol_u = solver_u.solve(b_u);
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            int k = (j - 1) * (nx - 1) + (i - 1);
            u_field[j][i] = sol_u(k);
        }
    }
}

void CavernaSolver2D::solve_poisson_v_A(VectorField& v_field, const VectorField& P) {
    int Nv = (nx) * (ny - 1);

    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx; ++i) {
            int k = (j - 1) * nx + (i - 1);
            b_v(k) = (P[j + 1][i] - P[j][i]) * hx; // Boundaries v=0 on all walls: no adjustment to b needed (implicit 0)
        }
    }
    sol_v = solver_v.solve(b_v);
    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx; ++i) {
            int k = (j - 1) * nx + (i - 1);
            v_field[j][i] = sol_v(k);
        }
    }
}