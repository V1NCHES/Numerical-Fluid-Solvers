#include "CavernaSolver2D.h"
#include "Utils.h"
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>

// ///////////////////////////////////////////////////////////////////////////
int CavernaSolver2D::run_full_algorithm(int max_iterations, double tolerance) 
{
    Timer timer;
    double t_sg;
    // u, v, p, tau, gamma, fx, fy, sigma = 0 По начальному приближению, а g = 0 из условий на границе
    //double diff_max = 1.0;
    for (; iter < max_iterations; ++iter) {
        
        // Step 1: Pressure Step (Empty block for your CG implementation)
        Timer timer_cg_custom;
        solve_pressure_cg_custom(); 
        t_sg = timer_cg_custom.elapsed();
        // Step 2: Update Gamma (Plasticity)
        update_gamma();
        
        // Step 3: Update Tau (Multiplier)
        diff_max.push_back(update_tau());
        
        // Step 4: Update F (Sources for next iter)
        update_f();

        std::cout << "Iteration " << iter << " | Diff Max: " << diff_max[iter- 1] << " | Time CG:" << t_sg << " s " << " | Time Algorithm: " << timer.elapsed() << "s" << std::endl;

        if (iter % 20 == 0)
        { 
            for (int j = 1; j <= ny - 1; ++j) {
                for (int i = 1; i <= nx - 1; ++i) {
                    f_psi[j][i] = (v[j][i + 1] - v[j][i]) / hx - (u[j + 1][i] - u[j][i]) / hy;
                }
            }
            solve_poisson_psi(psi, f_psi);
            save_iter(iter);
            save_all(iter);    
        }    
        
        if (diff_max[iter - 1] < tolerance) { std::cout << "Converged in " << iter<< " iterations." << std::endl; break;}
    }

    save_all(iter);
    save_iter(iter);
    return 0;
}

// ///////////////////////////////////////////////////////////////////////////
void CavernaSolver2D::norma_11_12(VectorField& q11, VectorField& q12, VectorField& q22)
{
    double T11, T22, T12;
    //norma 11
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            T11 = q11[j][i];  T22 = q22[j][i];  T12 = avg4(q12, i - 1, j - 1);

            norma11[j][i] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
             
            //ag_t12[j][i] = T12;
        }
    }

    for (int j = 1; j <= ny; ++j)
    {
        T11 = q11[j][0];   T22 = q22[j][0]; T12 = avg2y(q12, 0, j - 1);
        norma11[j][0] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
        //ag_t12[j][0] = T12;
         

        T11 = q11[j][nx + 1];   T22 = q22[j][nx + 1]; T12 = avg2y(q12, nx, j - 1);
        norma11[j][nx + 1] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
        //ag_t12[j][nx+1] = T12;
        
    }

    for (int i = 1; i <= nx; ++i)
    {
        T11 = q11[0][i];   T22 = q22[0][i]; T12 = avg2x(q12, i - 1, 0);
        norma11[0][i] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
        
        //ag_t12[0][i] = T12;
        T11 = q11[ny + 1][i];   T22 = q22[ny + 1][i]; T12 = avg2x(q12, i - 1, ny);
        norma11[ny + 1][i] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
        
        //ag_t12[ny+1][i] = T12;
    }


    // norma 12
    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            T11 = avg4(q11, i, j);  T22 = avg4(q22, i, j); T12 = q12[j][i];
            norma12[j][i] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
            

        }
    }

    for (int j = 1; j < ny; ++j)
    {
        T11 = avg2y(q11, 0, j);  T22 = avg2y(q22, 0, j); T12 = q12[j][0];
        norma12[j][0] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
        

        T11 = avg2y(q11, nx + 1, j);  T22 = avg2y(q22, nx + 1, j); T12 = q12[j][nx];
        norma12[j][nx] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
        
    }

    for (int i = 1; i < nx; ++i)
    {
        T11 = avg2x(q11, i, 0);  T22 = avg2x(q22, i, 0); T12 = q12[0][i];
        norma12[0][i] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
        

        T11 = avg2x(q11, i, ny + 1);  T22 = avg2x(q22, i, ny + 1); T12 = q12[ny][i];
        norma12[ny][i] = sqrt(T11 * T11 + T22 * T22 + 2 * T12 * T12);
        
    }

}

void CavernaSolver2D::update_gamma12(double modul, double dom, int i, int j)
{
    if (modul < tau_s)
        gamma12[j][i] = 0;
    else
        gamma12[j][i] = (1.0 - tau_s / modul) * t12[j][i] / dom;  
}

void CavernaSolver2D::update_gamma11_22(double modul, double dom, int i, int j)
{
    if (modul < tau_s)
    {
        gamma11[j][i] = 0; gamma22[j][i] = 0;
    }
    else
    {
        gamma11[j][i] = (1.0 - tau_s / modul) * t11[j][i] / dom; gamma22[j][i] = (1.0 - tau_s / modul) * t22[j][i] / dom;  
    }
}

void CavernaSolver2D::update_t()
{
    // gamma11, gamma22 depend on t11, t22 in centers
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i) {
            t11[j][i] = tau11[j][i] + 2.0 * r * du_dx[j][i];
            t22[j][i] = tau22[j][i] + 2.0 * r * dv_dy[j][i];
        }
    }
   // t11[0][0] = 0, t11[ny + 1][0] = 0, t11[0][nx + 1] = 0, t11[ny + 1][nx + 1] = 0;
    //t22[0][0] = 0, t22[ny + 1][0] = 0, t22[0][nx + 1] = 0, t22[ny + 1][nx + 1] = 0;

    // gamma12 depends on t12 in nodes
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) 
        {
            t12[j][i] = tau12[j][i] + r * (du_dy[j][i] + dv_dx[j][i]);
            D12[j][i] = (du_dy[j][i] + dv_dx[j][i])/2.0;
        }
    }
    //t12[0][0] = 0, t12[ny][0] = 0, t12[0][nx] = 0, t12[ny][nx] = 0;
    //D12[0][0] = 0, D12[ny][0] = 0, D12[0][nx] = 0, D12[ny][nx] = 0;
}

void CavernaSolver2D::update_gamma() 
{
    // Step 1: Вычисление производных dudx dvdy D12 = (dudy + dvdx)/2;
    compute_derivatives();
    double dom = 2.0 * (r + mu);
    double T11, T22, T12;
    // Step 2: Вычисление производных  t11 = tau11 + 2r dudx; t22 = tau22 + 2r dvdy;  t12 = tau12 + r D12; 
    update_t();
    // Step 3: Вычисление Нормы   norma11 для gamma11, gamma22 осреднение gamma12;  norma12 для gamma12 осреднение gamma11 gamma22
    norma_11_12(t11, t12, t22);

    // Step 4: Обновление gamma
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            update_gamma12(norma12[j][i], dom, i, j);
        }
    }

    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i) {
            update_gamma11_22(norma11[j][i], dom, i, j);
        }
    }
}

double CavernaSolver2D::update_tau() {

    double diff = 0.0;
    double T11, T22, T12;

    // Step 1: Обновление tau и вычисление diff 
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i) 
        {          
            T11 = (du_dx[j][i] - gamma11[j][i]);   tau11[j][i] += 2.0 * rho * T11;
            T22 = (dv_dy[j][i] - gamma22[j][i]);   tau22[j][i] += 2.0 * rho * T22;
                   
            diff = std::max({diff, std::abs(T11), std::abs(T22)});
        }
    }
    
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) 
        {
            T12 = (D12[j][i] - gamma12[j][i]);     tau12[j][i] += 2.0 * rho * T12;
            
            diff = std::max(diff, std::abs(T12));
        }
    }

    // Step 3: Вычисление Нормы   norma11 для tau11, tau22 осреднение tau12; norma12  для tau12 осреднение tau11 tau22
    norma_11_12(tau11, tau12, tau22);
    return diff;
}

void CavernaSolver2D::update_Q() {
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            Q11[j][i] = tau11[j][i] - 2.0 * r * gamma11[j][i];
            Q22[j][i] = tau22[j][i] - 2.0 * r * gamma22[j][i];
        }
    }
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            Q12[j][i] = tau12[j][i] - 2.0 * r * gamma12[j][i];
        }
    }
}

void CavernaSolver2D::update_f() 
{
    update_Q();

    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            fu[j][i] = (Q11[j][i+1] - Q11[j][i])/hx + (Q12[j][i] - Q12[j-1][i])/hy;
        }
    }
    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            fv[j][i] = (Q12[j][i] - Q12[j][i-1])/hx + (Q22[j+1][i] - Q22[j][i])/hy;
        }
    }
}

void CavernaSolver2D::compute_derivatives() {
    // 1. du_dx and dv_dy (Centers and ghost layers)
    
    
    // Internal centers
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            du_dx[j][i] = (u[j][i] - u[j][i - 1]) / hx;
            dv_dy[j][i] = (v[j][i] - v[j - 1][i]) / hy;
        }
    }
    double ux_0 = 0, ux_nx = 0;
    // du_dx boundary layers (i=0 and i=nx+1)
    for (int j = 1; j <= ny; ++j) {
        du_dx[j][0] = ( -3 * ux_0 + 4.0 * u[j][1] - u[j][2]) / (2.0 * hx); // Left wall u=0
        du_dx[j][nx + 1] = (3* ux_nx - 4.0 * u[j][nx - 1] + u[j][nx - 2]) / (2.0 * hx); // Right wall u=0
        //dv_dy[j][0] = 0 ; dv_dy[j][nx+1] =0 
    }
    double vy_0 = 0, vy_ny = 0;
    // dv_dy boundary layers (j=0 and j=ny+1)
    for (int i = 1; i <= nx; ++i) {
        dv_dy[0][i] = ( -3* vy_0 + 4.0 * v[1][i] - v[2][i]) / (2.0 * hy); // Bottom wall v=0
        dv_dy[ny + 1][i] = (3 * vy_ny - 4.0 * v[ny - 1][i] + v[ny - 2][i]) / (2.0 * hy); // Top wall v=0
        // du_dx[0][i] = 0 ; du_dx[ny + 1][i]] =0 
    }

    // 2. du_dy and dv_dx (Nodes)   
    // Internal nodes
    for (int j = 1; j <= ny - 1; ++j) {
        for (int i = 1; i <= nx - 1; ++i) {
            du_dy[j][i] = (u[j + 1][i] - u[j][i]) / hy;
            dv_dx[j][i] = (v[j][i + 1] - v[j][i]) / hx;
        }
    }

    double uy_0 = 0, uy_ny = 1.0;
    // Horizontal walls (j=0 and j=ny)
    for (int i = 1; i < nx; ++i) {
        //dv_dx[0][i] = 0.0;  dv_dx[ny][i] = 0.0;      
        du_dy[0][i] = (-8 * uy_0 + 9.0 * u[1][i] - u[2][i]) / (3.0 * hy); // Bottom wall u=0
        du_dy[ny][i] = (8 * uy_ny - 9.0 * u[ny][i] + u[ny - 1][i]) / (3.0 * hy); // Top wall u=1
    }
    double vx_0 = 0, vx_nx = 0;
    // Vertical walls (i=0 and i=nx)
    for (int j = 1; j < ny; ++j) {
        //du_dy[j][0] = 0.0;  du_dy[j][nx] = 0.0;
       
        dv_dx[j][0] = (-8 * vx_0 + 9.0 * v[j][1] - v[j][2]) / (3.0 * hx); // Left wall v=0
        dv_dx[j][nx] = (8 * vx_nx - 9.0 * v[j][nx] + v[j][nx - 1]) / (3.0 * hx); // Right wall v=0
    }
 
}

double CavernaSolver2D::avg4(const VectorField& field, int i, int j) {
    return 0.25 * (field[j][i] + field[j+1][i] + field[j][i+1] + field[j+1][i+1]);
}
double CavernaSolver2D::avg2x(const VectorField& field, int i, int j) {
    return 0.5 * (field[j][i] + field[j][i+1]);
}
double CavernaSolver2D::avg2y(const VectorField& field, int i, int j) {
    return 0.5 * (field[j][i] + field[j+1][i]);
}

