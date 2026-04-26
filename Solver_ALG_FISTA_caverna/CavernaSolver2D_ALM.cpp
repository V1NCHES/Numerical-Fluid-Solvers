#include"CavernaSolver2D.h"

int CavernaSolver2D::run_full_algorithm_boost(int max_iterations, double tolerance)
{
    Timer timer;
    double t_sg;
    //double diff_max = 1.0;
    // u, v, p, tau, gamma, fx, fy, sigma = 0 �� ���������� �����������, � g = 0 �� ������� �� �������
    if (iter == 1)
    //for (; iter < 41; ++iter)
    {
        //std::cout << "Iteration 0" << std::endl;
        // Step 0.1:
        Timer timer_cg_custom;
        solve_pressure_cg_custom();
        t_sg = timer_cg_custom.elapsed();
        // Step 0.2: Update Gamma (Plasticity)
        update_gamma();

        // Step 0.3: Update Tau (Multiplier)
        diff_max.push_back(update_tau());
        std::cout << "Iteration " << iter << " | Diff Max: " << diff_max[iter - 1] << " | Time CG:" << t_sg << " s " << " | Time Algorithm: " << timer.elapsed() << "s" << std::endl;
        iter++;
    }
    std::cout << std::endl << "Start fast ALM" << std::endl;
    for (; iter < max_iterations; ++iter) {

        // Step 1: Update Gamma (Plasticity)
        update_q();
        // Step 2: Update F (Sources for next iter)
        update_f_boost();
        // Step 3: Pressure Step (Empty block for your CG implementation)
        Timer timer_cg_custom;
        solve_pressure_cg_custom();
        t_sg = timer_cg_custom.elapsed();



        // Step 4: Update Tau (Multiplier)

        update_tau_boost();


        diff_max.push_back(extrapolate_tau_next_optimized()); // extrapolate_tau_next

        std::cout << "Iteration " << iter << " | Diff Max: " << diff_max[iter - 1] << " | Time CG:" << t_sg << " s " << " | Time Algorithm: " << timer.elapsed() << "s" << std::endl;

        //std::cout << "Iteration " << iter << std::endl << std::endl;
        if (iter % 10 == 0 ) //|| iter < 25
        {
            for (int j = 1; j <= ny - 1; ++j) {
               for (int i = 1; i <= nx - 1; ++i) {
                   f_psi[j][i] = (v[j][i + 1] - v[j][i]) / hx - (u[j + 1][i] - u[j][i]) / hy;
               }
            }
            solve_poisson_psi(psi, f_psi);

            save_all(iter);
            //save_iter_boost(iter);
        }
        if (diff_max[iter - 1] < tolerance) { std::cout << "Converged in " << iter << " iterations." << std::endl; break; }
    }

    
    
    save_all(iter);
   
    return 0;
}

double CavernaSolver2D::diff_l2_norma()
{
    double sum_sq = 0.0, diff11, diff12;
    int count = 0;

    // ��������� �������� ��� ���������� ���������� (11 � 22)
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i) {
            diff11 = norma11[j][i];
            sum_sq += diff11 * diff11;
            count++;
        }
    }

    // ��������� �������� ��� ����������� ���������� (12)
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            diff12 = norma12[j][i];
            sum_sq += diff12 * diff12;
            count++;
        }
    }

    // ���������� ������ �� �������� (��� Dl2 � MATLAB)
    return std::sqrt(sum_sq / count); //
}


double CavernaSolver2D::diff_max_norma()
{
    double diff = 0.0;
    double diff_11, diff_12;
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i)
        {
            diff_11 = norma11[j][i];
            diff = std::max({diff, std::abs(diff_11)});
        }
    }

    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i)
        {
            diff_12 = norma12[j][i];
            diff = std::max({ diff, std::abs(diff_12) });
        }
    }
    return diff;
}

void CavernaSolver2D::copy_tau(VectorField& Tau11, VectorField& Tau12, VectorField& Tau22, VectorField& f11, VectorField& f12, VectorField& f22)
{
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i)
        {
            f11[j][i] = Tau11[j][i];
            f22[j][i] = Tau22[j][i];
        }
    }

    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i)
        {
            f12[j][i] = Tau12[j][i];
        }
    }
}

double CavernaSolver2D::calculation_q(double tau_k, double norma_tau_k)
{
    if (norma_tau_k > tau_s)
        return (norma_tau_k - tau_s) * tau_k / norma_tau_k;
    return 0;
}



void CavernaSolver2D::update_q()
{
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i) {
            q11[j][i] = calculation_q(tau11[j][i], norma11[j][i]);
            q22[j][i] = calculation_q(tau22[j][i], norma11[j][i]);
        }
    }
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            q12[j][i] = calculation_q(tau12[j][i], norma12[j][i]);
        }
    }
}

void CavernaSolver2D::update_Q_boost()
{
    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            Q11[j][i] = tau11[j][i] -  q11[j][i];
            Q22[j][i] = tau22[j][i] -  q22[j][i];
        }
    }
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            Q12[j][i] = tau12[j][i] - q12[j][i];
        }
    }
}

void CavernaSolver2D::update_f_boost()
{
    update_Q_boost();

    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            fu[j][i] = (Q11[j][i + 1] - Q11[j][i]) / hx + (Q12[j][i] - Q12[j - 1][i]) / hy;
        }
    }
    for (int j = 1; j < ny; ++j) {
        for (int i = 1; i <= nx; ++i) {
            fv[j][i] = (Q12[j][i] - Q12[j][i - 1]) / hx + (Q22[j + 1][i] - Q22[j][i]) / hy;
        }
    }
}

void CavernaSolver2D::update_tau_boost()
{
    compute_derivatives();

    double diff = 0.0;
    double T11, T22, T12, tau11_prev, tau12_prev, tau22_prev;
    // Step 1: ���������� tau � ���������� diff 
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i)
        {
            tau11_prev = tau11_aux[j][i];       tau22_prev = tau22_aux[j][i];
            T11 = (du_dx[j][i] - q11[j][i]);   tau11_aux[j][i] = tau11[j][i] +  T11;
            T22 = (dv_dy[j][i] - q22[j][i]);   tau22_aux[j][i] = tau22[j][i] +  T22;

            diff = std::max({ diff, std::abs(T11), std::abs(T22) });
            diff_tau11[j][i] = tau11_aux[j][i] - tau11_prev; diff_tau22[j][i] = tau22_aux[j][i] - tau22_prev;
        }
    }

    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i)
        {
            tau12_prev = tau12_aux[j][i];
            D12[j][i] = (du_dy[j][i] + dv_dx[j][i]) / 2.0;
            T12 = (D12[j][i] - q12[j][i]);     tau12_aux[j][i] = tau12[j][i] +  T12;

            diff = std::max(diff, std::abs(T12));
            diff_tau12[j][i] = tau12_aux[j][i] - tau12_prev;
        }
    }
    norma_11_12(diff_tau11, diff_tau12, diff_tau22);
    //c_k = diff_max_norma();
    c_k = diff_l2_norma();
    //c_k = diff;
}




// Algorithm 9 (Fast AMA) + Algorithm 8 (Fast ADMM with Restart) -- Goldstein et al.
//
// Nesterov extrapolation on the dual variable (stress tensor):
//   alpha_k  = (1 + sqrt(1 + 4*alpha_{k-1}^2)) / 2       -- Alg.9 step 4
//   tau_new  = lambda_k + (alpha_{k-1}-1)/alpha_k * (lambda_k - lambda_{k-1})  -- Alg.9 step 5
//
//   lambda_k     = tau_aux       (current proximal step result)
//   lambda_{k-1} = tau_aux_prev  (previous proximal step result)
//
// Restart (Alg.8 step 11):  c_k > eta * c_{k-1}  OR  5-step monotone growth
// On restart: reset alpha=1, tau = lambda_k, lambda_{k-1} = lambda_k
double CavernaSolver2D::extrapolate_tau_next_optimized()
{
    double diff = 0.0;
    
    // eta_relax = 0.99 (��� � ������ Goldstein, Alg 8). 
    // ������� ������� �������� (0.999999) �������� "�������-�����" ��-�� ����.
    const double eta_relax = 0.9999999; 
    //double alpha_boost = 1.0, alpha_boost_prev = 1.0;
    // ���������� ����������� ������� ��� Cooldown (������ �������� ������� �����)
    
    double momentum = 0;
    // ������� ��������: ���� ������ � ������ ����� 5 ����� � �������� ��������
    if (c_k > eta_relax * c_k_prev)
    {
        std::cout << "; [Restart] c_k/c_prev=" << c_k / c_k_prev << " at iter " << iter << std::endl;
        
        
        alpha_boost_prev = 1.0; 
        c_k_prev = c_k / eta_relax; // ���������� ����� � �������� ��������

        for (int j = 0; j <= ny + 1; ++j) {
            for (int i = 0; i <= nx + 1; ++i) { 
                //tau11[j][i] = tau11_aux[j][i]; tau22[j][i] = tau22_aux[j][i];
                //tau11_aux_prev[j][i] = tau11_aux[j][i]; tau22_aux_prev[j][i] = tau22_aux[j][i];

                tau11[j][i] = tau11_aux_prev[j][i]; tau22[j][i] = tau22_aux_prev[j][i];
                tau11_aux[j][i] = tau11_aux_prev[j][i]; tau22_aux[j][i] = tau22_aux_prev[j][i];
            }
        }
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
               // tau12[j][i] = tau12_aux[j][i];  tau12_aux_prev[j][i] = tau12_aux[j][i]; 
                tau12_aux[j][i] = tau12_aux_prev[j][i];
                tau12[j][i] = tau12_aux_prev[j][i]; 
            }
        }
        //diff = c_k;
        diff = diff_max[iter-2];
    }
    else
    {
        alpha_boost = 0.5* (1.0 + sqrt(1.0 + 4.0 * alpha_boost_prev * alpha_boost_prev)) ;
        
       momentum = (alpha_boost_prev - 1.0) / alpha_boost;

        // �������������� ��� ���� (��� ���� ����������)
        for (int j = 0; j <= ny + 1; ++j) {
            for (int i = 0; i <= nx + 1; ++i) {
                double old11 = tau11[j][i], old22 = tau22[j][i];

                tau11[j][i] = tau11_aux[j][i] + momentum * (tau11_aux[j][i] - tau11_aux_prev[j][i]);  tau22[j][i] = tau22_aux[j][i] + momentum * (tau22_aux[j][i] - tau22_aux_prev[j][i]);

                tau11_aux_prev[j][i] = tau11_aux[j][i]; tau22_aux_prev[j][i] = tau22_aux[j][i];

                diff_tau11[j][i] = tau11[j][i] - old11;  diff_tau22[j][i] = tau22[j][i] - old22;
               
            }
        }

        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                double old12 = tau12[j][i];
                tau12[j][i] = tau12_aux[j][i] + momentum * (tau12_aux[j][i] - tau12_aux_prev[j][i]);
                tau12_aux_prev[j][i] = tau12_aux[j][i];
                diff_tau12[j][i] = tau12[j][i] - old12;
            }
        }

        alpha_boost_prev = alpha_boost;
        c_k_prev = c_k;

        norma_11_12(diff_tau11, diff_tau12, diff_tau22);
        //diff = diff_l2_norma();
        diff = diff_max_norma();
    }


    norma_11_12(tau11, tau12, tau22);
    return diff;
}

