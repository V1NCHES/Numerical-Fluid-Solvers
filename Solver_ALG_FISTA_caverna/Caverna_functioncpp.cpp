#include "CavernaSolver2D.h"

CavernaSolver2D::CavernaSolver2D(double h, double Lx, double Ly, double mu, double tau_s, double r, double rho, const std::string& output_dir)
    : hx(h), hy(h), Lx(Lx), Ly(Ly), mu(mu), tau_s(tau_s), r(r), rho(rho), output_dir(output_dir) {

    nx = static_cast<int>(Lx / h);
    ny = static_cast<int>(Ly / h);
    //h2 = hx * hx;
    std::cout << "nx = " << nx << " ; ny = " << ny << std::endl;
    h2 = hx * hx;
    // Initialize thread pool with 2 threads for U and V solvers
    pool = std::make_unique<ThreadPool>(2);

    init_fields();
    precompute_matrices();
}

void CavernaSolver2D::boundary_conditions()
{
    double Pi = acos(-1);
    double K = 1.0 / (4 * Pi * Pi), k = 1.0 / Pi, pi2 = 2 * Pi, pi_2 = Pi * Pi;
    /*for (int i = 0; i <= nx_d; ++i)
        v1[i][ny1] = -0.1;*/
    double x = 0.0, y = 0.0;
    double hx_2 = hx / 2.0, hy_2 = hy / 2.0;
    // std::cout << "hx = " << hx << " ; hy = " << hy << "hx_2 = " << hx_2 << " ; hy_2 = " << hy_2 << "\n";
    for (int j = 1; j < ny + 1; ++j)
        for (int i = 1; i < nx; ++i)
        {
            x = i * hx, y = (j)*hy - hy_2;

            // std::cout << "x = " << x << " ; y = " << y << "\n";
            fu[j][i] = sin(pi2 * y);
            //u1[j][i] = K * (1 - cos(pi2 * x)) * sin(pi2 * y);

        }

    // std::cout << "hx = " << hx << " ; hy = " << hy << "hx_2 = " << hx_2 << " ; hy_2 = " << hy_2 << "\n";
    for (int j = 1; j < ny; ++j)
        for (int i = 1; i < nx + 1; ++i)
        {
            x = (i)*hx - hx_2, y = j * hy;
            // std::cout << "x = " << x << " ; y = " << y << "\n";
            fv[j][i] = sin(pi2 * x);
            //v1[j][i] = K * sin(pi2 * x) * (1 - cos(pi2 * y));

        }
    //std::cout << "hx = " << hx << " ; hy = " << hy << "hx_2 = " << hx_2 << " ; hy_2 = " << hy_2 << "\n";
    for (int j = 1; j < ny + 1; ++j)
        for (int i = 1; i < nx + 1; ++i)
        {
            x = (i)*hx - hx_2, y = (j)*hy - hy_2;
            // std::cout << "x = " << x << " ; y = " << y << "\n";
            g[j][i] = -k * sin(pi2 * x) * sin(pi2 * y);
        }
}

void CavernaSolver2D::symmetric_x(const VectorField& field, int nx, int ny)
{
    double eps = 1e-10;
    bool is_symmetric = true;

    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx / 2; ++i) {
            double val1 = field[j][i];
            double val2 = field[j][nx - i];
            if (std::abs(val1 - val2) > eps) {
                std::cout << "[Symmetry X FAIL] at j=" << j << ", i=" << i
                    << " (left: " << val1 << ", right: " << val2
                    << ", diff: " << std::abs(val1 - val2) << ")" << std::endl;
                is_symmetric = false;
                return; // Выходим при первом нарушении
            }
        }
    }

    if (is_symmetric) {
        std::cout << "\n[Symmetry X OK] Field is symmetric relative to X axis." << std::endl;
    }
}

void CavernaSolver2D::antysymetric_x(const VectorField& field, int nx, int ny)
{
    double eps = 1e-10;
    bool is_antisymmetric = true;

    for (int j = 1; j <= ny; ++j) {
        for (int i = 1; i <= nx / 2; ++i) {
            double val1 = field[j][i];
            double val2 = field[j][nx - i];
            if (std::abs(val1 + val2) > eps) {
                std::cout << "[Antisymmetry X FAIL] at j=" << j << ", i=" << i
                    << " (left: " << val1 << ", right: " << val2
                    << ", sum: " << std::abs(val1 + val2) << ")" << std::endl;
                is_antisymmetric = false;
                return; // Выходим при первом нарушении
            }
        }
    }

    if (is_antisymmetric) {
        std::cout << "\n[Antisymmetry X OK] Field is antisymmetric relative to X axis." << std::endl;
    }
}


void CavernaSolver2D::init_fields() {
    int rows = ny + 2;
    int cols = nx + 2;

    u.assign(rows, std::vector<double>(cols, 0.0));
    v.assign(rows, std::vector<double>(cols, 0.0));
    p.assign(rows, std::vector<double>(cols, 0.0));
    psi.assign(rows, std::vector<double>(cols, 0.0));

    u_tilde.assign(rows, std::vector<double>(cols, 0.0));
    v_tilde.assign(rows, std::vector<double>(cols, 0.0));
    r1.assign(rows, std::vector<double>(cols, 0.0));
    Ap.assign(rows, std::vector<double>(cols, 0.0));
    //Bx_p.assign(rows, std::vector<double>(cols, 0.0));
    //By_p.assign(rows, std::vector<double>(cols, 0.0));

    p_prev.assign(rows, std::vector<double>(cols, 0.0));
    F.assign(rows, std::vector<double>(cols, 0.0));
    g.assign(rows, std::vector<double>(cols, 0.0));
    //BxT.assign(rows, std::vector<double>(cols, 0.0));
    //ByT.assign(rows, std::vector<double>(cols, 0.0));

    //d1.assign(rows, std::vector<double>(cols, 0.0));

    tau11.assign(rows, std::vector<double>(cols, 0.0));
    tau12.assign(rows, std::vector<double>(cols, 0.0));
    tau22.assign(rows, std::vector<double>(cols, 0.0));

    gamma11.assign(rows, std::vector<double>(cols, 0.0));
    gamma12.assign(rows, std::vector<double>(cols, 0.0));
    gamma22.assign(rows, std::vector<double>(cols, 0.0));

    Q11.assign(rows, std::vector<double>(cols, 0.0));
    Q12.assign(rows, std::vector<double>(cols, 0.0));
    Q22.assign(rows, std::vector<double>(cols, 0.0));


    norma12.assign(rows, std::vector<double>(cols, 0.0));
    norma11.assign(rows, std::vector<double>(cols, 0.0));

    fu.assign(rows, std::vector<double>(cols, 0.0));
    fv.assign(rows, std::vector<double>(cols, 0.0));
    f_psi.assign(rows, std::vector<double>(cols, 0.0));


    du_dx.assign(rows, std::vector<double>(cols, 0.0));
    du_dy.assign(rows, std::vector<double>(cols, 0.0));
    dv_dx.assign(rows, std::vector<double>(cols, 0.0));
    dv_dy.assign(rows, std::vector<double>(cols, 0.0));
    D12.assign(rows, std::vector<double>(cols, 0.0));
    

    t11.assign(rows, std::vector<double>(cols, 0.0));
    t12.assign(rows, std::vector<double>(cols, 0.0));
    t22.assign(rows, std::vector<double>(cols, 0.0));

    //BxT.assign(rows, std::vector<double>(cols, 0.0));
    //ByT.assign(rows, std::vector<double>(cols, 0.0));


    q11.assign(rows, std::vector<double>(cols, 0.0));
    q12.assign(rows, std::vector<double>(cols, 0.0));
    q22.assign(rows, std::vector<double>(cols, 0.0));


    tau11_aux_prev.assign(rows, std::vector<double>(cols, 0.0));
    tau12_aux_prev.assign(rows, std::vector<double>(cols, 0.0));
    tau22_aux_prev.assign(rows, std::vector<double>(cols, 0.0));

    tau11_aux.assign(rows, std::vector<double>(cols, 0.0));
    tau12_aux.assign(rows, std::vector<double>(cols, 0.0));
    tau22_aux.assign(rows, std::vector<double>(cols, 0.0));

    tau11_aux_stable.assign(rows, std::vector<double>(cols, 0.0));
    tau12_aux_stable.assign(rows, std::vector<double>(cols, 0.0));
    tau22_aux_stable.assign(rows, std::vector<double>(cols, 0.0));

    

    diff_tau11.assign(rows, std::vector<double>(cols, 0.0));
    diff_tau12.assign(rows, std::vector<double>(cols, 0.0));
    diff_tau22.assign(rows, std::vector<double>(cols, 0.0));
    
    //ag_t12.assign(rows, std::vector<double>(cols, 0.0));

    //eta = 0.999;
    eta = 1.0;
    alpha_boost_prev = 1.0;
    c_k_prev = 1000;
    // ////////////////////////////
    iter = 1;
    next_boost = false;

    //update_iter(iter); iter++;
}


void CavernaSolver2D::save_all(int iteration) {
    std::string prefix = output_dir + "/iret" + std::to_string(iteration) + "_";
    //std::string prefix = output_dir + "/";

    //if (iteration == 1)

   /* save_field(prefix + "u_iteration_1.txt", u_tilde, nx, ny, hx, hy, 1);
        save_field(prefix + "v_iteration_1.txt", v_tilde, nx, ny, hx, hy, 2);*/
        //save_field(prefix + "p_iteration_1.txt", p, nx, ny, hx, hy, 3);

        //save_field(prefix + "psi.txt", psi, nx, ny, hx, hy, 4);
         //save_field(prefix + "u.txt", u, nx, ny, hx, hy, 1);
         //save_field(prefix + "v.txt", v, nx, ny, hx, hy, 2);
         //save_field(prefix + "p.txt", p, nx, ny, hx, hy, 3);


         ////save_field(prefix + "u.txt", u, nx, ny, hx, hy, 1);
         ////save_field(prefix + "v.txt", v, nx, ny, hx, hy, 2);
         ////save_field(prefix + "p.txt", p, nx, ny, hx, hy, 3);
    save_field(prefix + "psi.txt", psi, nx, ny, hx, hy, 4);

    //// Save stresses (tau)
    /*save_field(prefix + "tau11.txt", tau11, nx, ny, hx, hy, 5);
    save_field(prefix + "tau22.txt", tau22, nx, ny, hx, hy, 5);
    save_field(prefix + "tau12.txt", tau12, nx, ny, hx, hy, 4);*/

    ////// Save deformation rates (gamma)
    //save_field(prefix + "gamma11.txt", gamma11, nx, ny, hx, hy, 5);
    //save_field(prefix + "gamma22.txt", gamma22, nx, ny, hx, hy, 5);
    //save_field(prefix + "gamma12.txt", gamma12, nx, ny, hx, hy, 4);

    save_field(prefix + "norma12.txt", norma12, nx, ny, hx, hy, 4);
    save_field(prefix + "norma11.txt", norma11, nx, ny, hx, hy, 5);
   
    /*save_field(prefix + "dudx.txt", du_dx, nx, ny, hx, hy, 5);
    save_field(prefix + "dvdy.txt", dv_dy, nx, ny, hx, hy, 5);
    save_field(prefix + "D12.txt", D12, nx, ny, hx, hy, 4);*/

    std::string prefix_1 = output_dir + "/";
    std::string filename = prefix_1 + "diff_max.txt";
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for saving: " << filename << std::endl;
        return ;
    }
    file << std::fixed << std::setprecision(8);

    for (int i = 0; i < diff_max.size(); ++i)
        file << diff_max[i] << std::endl;
    file.close();
}



void CavernaSolver2D::save_iter_boost(int iteration)
{
    std::string prefix = "SAVE ITERATION Boost/iter_" + std::to_string(iteration) + "_";
    int NX = nx + 2, Ny = ny + 2;
    //save_iteration(prefix + "fu.txt", fu, NX, Ny);
    //save_iteration(prefix + "fv.txt", fv, NX, Ny);

    save_iteration(prefix + "tau11.txt", tau11, NX, Ny);
    save_iteration(prefix + "tau22.txt", tau22, NX, Ny);
    save_iteration(prefix + "tau12.txt", tau12, NX, Ny);

    save_iteration(prefix + "tau11_aux_prev.txt", tau11_aux_prev, NX, Ny);
    save_iteration(prefix + "tau12_aux_prev.txt", tau12_aux_prev, NX, Ny);
    save_iteration(prefix + "tau22_aux_prev.txt", tau22_aux_prev, NX, Ny);

    

    save_iteration(prefix + "tau11_aux.txt", tau11_aux, NX, Ny);
    save_iteration(prefix + "tau12_aux.txt", tau12_aux, NX, Ny);
    save_iteration(prefix + "tau22_aux.txt", tau22_aux, NX, Ny);

    save_iteration(prefix + "q11.txt", q11, NX, Ny);
    save_iteration(prefix + "q22.txt", q22, NX, Ny);
    save_iteration(prefix + "q12.txt", q12, NX, Ny);

    save_iteration(prefix + "norma11.txt", norma11, NX, Ny);
    save_iteration(prefix + "norma12.txt", norma12, NX, Ny);

    std::ofstream file("SAVE ITERATION/c_k_diff_Max.txt");
    if (!file.is_open()) {
        std::cerr << "Error opening file for saving: " << "SAVE ITERATION/c_k_diff_Max.txt" << std::endl;
        return;
    }
    file << std::fixed << std::setprecision(8);

    file << c_k_prev << " " << c_k << std::endl;

    for (size_t j = 0; j < diff_max.size(); ++j)
    {
        file << j << " " << diff_max[j] << std::endl;
    }
    file.close();
}

void CavernaSolver2D::update_iter_boost(int iteration)
{
    std::string prefix = "SAVE ITERATION/iter_" + std::to_string(iteration) + "_";

    //reading(prefix + "fu.txt", fu);
    //reading(prefix + "fv.txt", fv);

    reading(prefix + "tau11.txt", tau11);
    reading(prefix + "tau22.txt", tau22);
    reading(prefix + "tau12.txt", tau12);

    reading(prefix + "tau11_aux_prev.txt", tau11_aux_prev);
    reading(prefix + "tau12_aux_prev.txt", tau12_aux_prev);
    reading(prefix + "tau22_aux_prev.txt", tau22_aux_prev);

    
    reading(prefix + "tau11_aux.txt", tau11_aux);
    reading(prefix + "tau12_aux.txt", tau12_aux);
    reading(prefix + "tau22_aux.txt", tau22_aux);

    reading(prefix + "q11.txt", q11);
    reading(prefix + "q22.txt", q22);
    reading(prefix + "q12.txt", q12);

    reading(prefix + "norma11.txt", norma11);
    reading(prefix + "norma12.txt", norma12);
}








void CavernaSolver2D::save_iter(int iteration)
{
    std::string prefix = "SAVE ITERATION/iter_" + std::to_string(iteration) + "_";
    int NX = nx + 2, Ny = ny + 2;
    save_iteration(prefix + "fu.txt", fu, NX, Ny);
    save_iteration(prefix + "fv.txt", fv, NX, Ny);

    save_iteration(prefix + "tau11.txt", tau11, NX, Ny);
    save_iteration(prefix + "tau22.txt", tau22, NX, Ny);
    save_iteration(prefix + "tau12.txt", tau12, NX, Ny);

    save_iteration(prefix + "gamma11.txt", gamma11, NX, Ny);
    save_iteration(prefix + "gamma22.txt", gamma22, NX, Ny);
    save_iteration(prefix + "gamma12.txt", gamma12, NX, Ny);
}

void CavernaSolver2D::update_iter(int iteration)
{
    std::string prefix = "SAVE ITERATION/iter_" + std::to_string(iteration) + "_";

    reading(prefix + "fu.txt", fu);
    reading(prefix + "fv.txt", fv);

    reading(prefix + "tau11.txt", tau11);
    reading(prefix + "tau22.txt", tau22);
    reading(prefix + "tau12.txt", tau12);

    reading(prefix + "gamma11.txt", gamma11);
    reading(prefix + "gamma22.txt", gamma22);
    reading(prefix + "gamma12.txt", gamma12);
}
