#include"CavernaSolver2D.h"

int CavernaSolver2D::run_full_algorithm_boost(int max_iterations, double tolerance)
{
    Timer timer;
    double t_sg;
    //double diff_max = 1.0;
    // u, v, p, tau, gamma, fx, fy, sigma = 0 По начальному приближению, а g = 0 из условий на границе
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
       // c_k_prev = diff_max[iter - 1];
        //update_f();

        std::cout << "Iteration " << iter << " | Diff Max: " << diff_max[iter - 1] << " | Time CG:" << t_sg << " s " << " | Time Algorithm: " << timer.elapsed() << "s" << std::endl;

        /*if (iter % 20 == 0)
        {
            for (int j = 1; j <= ny - 1; ++j) {
                for (int i = 1; i <= nx - 1; ++i) {
                    f_psi[j][i] = (v[j][i + 1] - v[j][i]) / hx - (u[j + 1][i] - u[j][i]) / hy;
                }
            }
            solve_poisson_psi(psi, f_psi);
            save_iter(iter);
            save_all(iter);
        }*/

        //std::cout << "Iteration " << iter << " | Diff Max: " << diff_max[iter-1] << " | Time CG:" << t_sg << " s " << " | Time Algorithm: " << timer.elapsed() << "s" << std::endl;
        //save_all(1);
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
        //if(diff_max[iter - 1])
        //if (diff_max[iter - 1] < tolerance) { std::cout << "Converged in " << iter << " iterations." << std::endl; break; }
    }

    
    
    save_all(iter);
    //save_iter_boost(iter);
    return 0;
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
    // Step 1: Обновление tau и вычисление diff 
    for (int j = 0; j <= ny + 1; ++j) {
        for (int i = 0; i <= nx + 1; ++i)
        {
            tau11_prev = tau11_aux[j][i];       tau22_prev = tau22_aux[j][i];
            T11 = (du_dx[j][i] - q11[j][i]);   tau11_aux[j][i] = tau11[j][i] + r * T11;
            T22 = (dv_dy[j][i] - q22[j][i]);   tau22_aux[j][i] = tau22[j][i] + r * T22;

            diff = std::max({ diff, std::abs(T11), std::abs(T22) });
            diff_tau11[j][i] = tau11_aux[j][i] - tau11_prev; diff_tau22[j][i] = tau22_aux[j][i] - tau22_prev;
        }
    }

    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i)
        {
            tau12_prev = tau12_aux[j][i];
            D12[j][i] = (du_dy[j][i] + dv_dx[j][i]) / 2.0;
            T12 = (D12[j][i] - q12[j][i]);     tau12_aux[j][i] = tau12[j][i] + r * T12;

            diff = std::max(diff, std::abs(T12));
            diff_tau12[j][i] = tau12_aux[j][i] - tau12_prev;
        }
    }
    norma_11_12(diff_tau11, diff_tau12, diff_tau22);
    c_k = diff_max_norma();
    //c_k = diff;
}

double CavernaSolver2D::extrapolate_tau_next()
{
    // double tau = 0;
    double diff = 0.0;
    double tau_11_cur, tau_12_cur, tau_22_cur;
    std::cout << std::endl<< "c_k = " << c_k << " ; eta * c_k_prev = " << eta * c_k_prev;


    if (c_k < eta * c_k_prev ) //|| next_boost
    {
        next_boost = false;
        alpha_boost = (1 + sqrt(1 + 4 * alpha_boost_prev * alpha_boost_prev)) / 2.0;
        std::cout << ";  boost(); alpha_boost = " << alpha_boost << std::endl;
        for (int j = 0; j <= ny + 1; ++j) {
            for (int i = 0; i <= nx + 1; ++i)
            {
                tau_11_cur = tau11[j][i]; tau_22_cur = tau22[j][i];
                tau11[j][i] = tau11_aux[j][i] + (alpha_boost_prev - 1.0) * (tau11_aux[j][i] - tau11_aux_prev[j][i]) / alpha_boost;

                tau22[j][i] = tau22_aux[j][i] + (alpha_boost_prev - 1.0) * (tau22_aux[j][i] - tau22_aux_prev[j][i]) / alpha_boost;

                tau11_aux_stable[j][i] = tau11_aux_prev[j][i];  tau11_aux_prev[j][i] = tau11_aux[j][i];

                tau22_aux_stable[j][i] = tau22_aux_prev[j][i]; tau22_aux_prev[j][i] = tau22_aux[j][i];

                //tau11_boost_prev[j][i] = tau11_boost[j][i];
                //tau22_boost_prev[j][i] = tau22_boost[j][i];
                diff_tau11[j][i] = tau11[j][i] - tau_11_cur;
                diff_tau22[j][i] = tau22[j][i] - tau_22_cur;
            }
        }

        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i)
            {
                tau_12_cur = tau12[j][i];

                tau12[j][i] = tau12_aux[j][i] + (alpha_boost_prev - 1.0) * (tau12_aux[j][i] - tau12_aux_prev[j][i]) / alpha_boost;

                tau12_aux_stable[j][i] = tau12_aux_prev[j][i];

                tau12_aux_prev[j][i] = tau12_aux[j][i];
                 //tau12_boost_prev[j][i] = tau12_boost[j][i];
                diff_tau12[j][i] = tau12[j][i] - tau_12_cur;
            }
        }

        alpha_boost_prev = alpha_boost;
        c_k_prev = c_k;
        norma_11_12(diff_tau11, diff_tau12, diff_tau22);
        diff = diff_max_norma();
    }
    else
    {
        if (iter > 3)
        {
            //std::cout << "()" << std::endl;
            //if (next_boost)
            //{
            //    //next_boost = false;
            //    alpha_boost_prev = 1.1;
            //}
            //else
            //{
            //    alpha_boost_prev = 1.0;
            //    next_boost = true;
            //}
            alpha_boost_prev = 1.0; //(1 + sqrt(5.0)) / 2.0;;
            next_boost = true;
            c_k_prev = c_k_prev / eta;
            std::cout << "; Anty_boost(); alpha_boost_prev = " << alpha_boost_prev << "; diff = "<< diff_max[iter - 2] << std::endl;

            for (int j = 0; j <= ny + 1; ++j) {
                for (int i = 0; i <= nx + 1; ++i)
                {
                    tau11[j][i] = tau11_aux_stable[j][i];

                    tau22[j][i] = tau22_aux_stable[j][i];

                    tau11_aux_prev[j][i] = tau11_aux_stable[j][i];
                    tau22_aux_prev[j][i] = tau22_aux_stable[j][i];
                    //tau11_boost[j][i] = tau11_boost_prev[j][i];

                    //tau22_boost[j][i] = tau22_boost_prev[j][i];
                }
            }

            for (int j = 0; j <= ny; ++j) {
                for (int i = 0; i <= nx; ++i)
                {

                    tau12[j][i] = tau12_aux_stable[j][i];
                    tau12_aux_prev[j][i] = tau12_aux_stable[j][i];
                    //tau12_boost[j][i] = tau12_boost_prev[j][i];
                }
            }
           
            diff = diff_max[iter - 2];
            
        }
        else
        {
            std::cout << "Error()" << std::endl;
            //return update_tau();

        }
        
    }

    norma_11_12(tau11, tau12, tau22);

    return diff;
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
    
    // eta_relax = 0.99 (как в статье Goldstein, Alg 8). 
    // Слишком большие значения (0.999999) вызывают "рестарт-петлю" из-за шума.
    const double eta_relax = 0.99; 

    // Используем статический счетчик для Cooldown (запрет рестарта слишком часто)
    static int last_restart_iter = -10; 
    double momentum = 0;
    // Условие рестарта: рост ошибки И прошло более 5 шагов с прошлого рестарта
    if (c_k > eta_relax * c_k_prev && (iter - last_restart_iter) > 5)
    {
        std::cout << "; [Restart] c_k/c_prev=" << c_k / c_k_prev << " at iter " << iter << std::endl;
        last_restart_iter = iter;
        
        alpha_boost_prev = 1.0; 
        c_k_prev = c_k; // Сбрасываем порог к текущему значению

        for (int j = 0; j <= ny + 1; ++j) {
            for (int i = 0; i <= nx + 1; ++i) {
                tau11[j][i] = tau11_aux[j][i];
                tau22[j][i] = tau22_aux[j][i];
                tau11_aux_prev[j][i] = tau11_aux[j][i];
                tau22_aux_prev[j][i] = tau22_aux[j][i];
            }
        }
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                tau12[j][i] = tau12_aux[j][i];
                tau12_aux_prev[j][i] = tau12_aux[j][i];
            }
        }
        diff = c_k;
    }
    else
    {
        alpha_boost = (1.0 + sqrt(1.0 + 4.0 * alpha_boost_prev * alpha_boost_prev)) / 2.0;
        if(iter > 50)
            momentum = 0.8* (alpha_boost_prev - 1.0) / alpha_boost;
        else
            momentum = (alpha_boost_prev - 1.0) / alpha_boost;

        // Экстраполируем ВСЕ узлы (как было изначально)
        for (int j = 0; j <= ny + 1; ++j) {
            for (int i = 0; i <= nx + 1; ++i) {
                double old11 = tau11[j][i];
                double old22 = tau22[j][i];

                tau11[j][i] = tau11_aux[j][i] + momentum * (tau11_aux[j][i] - tau11_aux_prev[j][i]);
                tau22[j][i] = tau22_aux[j][i] + momentum * (tau22_aux[j][i] - tau22_aux_prev[j][i]);

                tau11_aux_prev[j][i] = tau11_aux[j][i];
                tau22_aux_prev[j][i] = tau22_aux[j][i];

                diff_tau11[j][i] = tau11[j][i] - old11;
                diff_tau22[j][i] = tau22[j][i] - old22;
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
        diff = diff_max_norma();
    }

    // Это ГЛАВНОЕ для фикса артефактов в углах: 
    // Причесываем границы после математической экстраполяции
    boundary_conditions(); 

    norma_11_12(tau11, tau12, tau22);
    return diff;
}

