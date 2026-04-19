#ifndef CAVERNA_SOLVER_2D_H
#define CAVERNA_SOLVER_2D_H

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <memory>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <omp.h>
#include "Utils.h"

using VectorField = std::vector<std::vector<double>>;

class CavernaSolver2D {
public:
    CavernaSolver2D(double h, double Lx, double Ly, double mu, double tau_s, double r, double rho, const std::string& output_dir);

    // Main Algorithm Loop
    int run_full_algorithm(int max_iterations, double tolerance);

    // Save results to file
    void save_all(int iteration);
    void save_iter(int iteration);
    void update_iter(int iteration);

    // Poisson Solvers (Public for user interaction if needed)
    void solve_poisson_u(VectorField& u_field, const VectorField& rhs_field);
    void solve_poisson_v(VectorField& v_field, const VectorField& rhs_field);
    void solve_poisson_psi(VectorField& psi_field, const VectorField& rhs_field);

    void solve_poisson_u_p(VectorField& u_field, const VectorField& rhs_field, VectorField& p_field);
    void solve_poisson_v_p(VectorField& v_field, const VectorField& rhs_field, VectorField& p_field);

    // Pressure Operator for CG
   
    void solve_pressure_cg_custom();
    void Operator_A(const VectorField& p_field);
    double scalar_product();
    void calculation_r();
    void calculation_g();

//public:
    // Parameters
    double hx, hy, h2;
    double Lx, Ly;
    int nx, ny;
    int iter;
    double mu, tau_s, r, rho;
    double r_cur, r_prev, tau_cur, tau_prev, alpha_cur, alpha_prev;
    //double diff_max;
    std::string output_dir;

    // Main Fields
    VectorField u, v, p, psi;
    VectorField tau11, tau12, tau22;
    VectorField norma12, norma11;
    VectorField gamma11, gamma12, gamma22;
    VectorField Q11, Q12, Q22;
    VectorField fu, fv, f_psi; // Sources (divergence of stresses)
    //void test();

    // Auxiliary Fields for CG and intermediate steps
    VectorField u_tilde, v_tilde, r1, Ap, Bx_p, By_p;
    VectorField p_prev, F, g; // CG vectors

    // Intermediate fields for derivatives
    VectorField du_dx, du_dy, dv_dx, dv_dy, D12;
    VectorField t11, t12, t22;

    // Eigen solver members
    Eigen::SparseMatrix<double> Au, Av, Apsi, Bxu, BxTu, Byv, ByTv;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver_u, solver_v, solver_psi;
    Eigen::VectorXd b_u, b_v, b_psi;
    Eigen::VectorXd sol_u, sol_v, sol_psi;

    // Thread pool for parallel execution
    std::unique_ptr<ThreadPool> pool;

    // Initialization
    void init_fields();
    void precompute_matrices();

    // Tensor/Algorithm steps
    void update_gamma();
    void update_t();
    void update_gamma12(double modul, double dom, int i, int j);
    void update_gamma11_22(double modul, double dom, int i, int j);

    double update_tau();
    void norma_11_12(VectorField& q11, VectorField& q12, VectorField& q22);

    void update_Q();
    void update_f();

    // Derivatives
    void compute_derivatives();

    // Averaging helpers
    double avg4(const VectorField& field, int i, int j);
    double avg2x(const VectorField& field, int i, int j);
    double avg2y(const VectorField& field, int i, int j);

    // symmetric x, antysymetric x
    void symmetric_x(const VectorField& field, int nx, int ny);
    void antysymetric_x(const VectorField& field, int nx, int ny);

    void boundary_conditions();
    void solve_poisson_u_caverna(VectorField& u_field, const VectorField& rhs_field);

    void solve_poisson_u_A(VectorField& u_field, const VectorField& rhs_field);
    void solve_poisson_v_A(VectorField& v_field, const VectorField& rhs_field);


    /// ////////////////////
    void copy_tau(VectorField& Tau11, VectorField& Tau12, VectorField& Tau22, VectorField& f11, VectorField& f12, VectorField& f22);

    VectorField q11, q12, q22;
    void update_q();
    double calculation_q(double tau_k, double norma_tau_k);
    void update_Q_boost();
    void update_f_boost();

    void update_tau_boost();
    VectorField  tau11_aux_prev, tau12_aux_prev, tau22_aux_prev, tau11_aux, tau12_aux, tau22_aux, tau11_aux_stable, tau12_aux_stable, tau22_aux_stable;

    double c_k, c_k_prev, eta;
    double alpha_boost_prev, alpha_boost;
    double extrapolate_tau_next();
    int run_full_algorithm_boost(int max_iterations, double tolerance);
    std::vector<double> diff_max;
    void save_iter_boost(int iteration);
    void update_iter_boost(int iteration);
    bool next_boost;

    /////////////////////////////////////

    double diff_max_norma();
    VectorField diff_tau11, diff_tau12, diff_tau22;
    double extrapolate_tau_next_optimized();
    //VectorField ag_t11, ag_t12;
};


#endif
