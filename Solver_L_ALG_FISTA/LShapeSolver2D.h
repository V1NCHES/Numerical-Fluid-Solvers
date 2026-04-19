#pragma once
#include <vector>
#include <iostream>
#include <map>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <omp.h>           // for parallel for loops (non-task parallelism)
#include <algorithm>       // for std::max
#include "ThreadPool.h"    // pre-allocated thread pool for tasks

typedef std::vector<std::vector<double>> VectorField;

struct LShapeConfig {
    int nx1, nx2, ny1, ny2;
    double h;
    double v_inflow;
    int lx_hole;        // Number of x-cells where top inflow acts
    std::string output_dir = "output_data_boost";  // Output directory
};

class LShapeSolver2D {
public:
    LShapeSolver2D(LShapeConfig config);
    ~LShapeSolver2D() = default;

    // Bingham parameters
    double mu = 1.0;
    double r_param = 1.0; // r in uzawa algorithm
    double rho = 1.0;
    double tau_s = 0.0;   // To be set from main()
    
    int iter = 0;
    std::vector<double> diff_max;

    int run_full_algorithm(int max_iterations, double tolerance);

    // Bingham algorithm steps
    void update_gamma();
    double update_tau();
    void update_f();
    void compute_derivatives();
    void update_t();
    void update_gamma12(double modul, double dom, int i, int j);
    void update_gamma11_22(double modul, double dom, int i, int j);
    double avg4(const VectorField& field, int j, int i);
    double avg2x(const VectorField& field, int j, int i);
    double avg2y(const VectorField& field, int j, int i);
    void norma_11_12(const VectorField& q11, const VectorField& q12, const VectorField& q22);


    // Preparation
    void prepare_fields();
    void precompute_matrices();

    // Core solvers
    void solve_poisson_u(VectorField& u_field, const VectorField& rhs_field);
    void solve_poisson_v(VectorField& v_field, const VectorField& rhs_field);

    // Operator A: solves Laplace for u_tilde (rhs = Bx*p) and v_tilde (rhs = By*p)
    //             then computes Ap = -(div u_tilde)
    void Operator_A(const VectorField& p_field);

    // Helpers for Operator_A  (Neumann ghost for pressure)
    void solve_poisson_u_A(VectorField& u_field, const VectorField& p_field);
    void solve_poisson_v_A(VectorField& v_field, const VectorField& p_field);

    // Final velocity correction:  -Delta u = fu - Bx*p,  -Delta v = fv - By*p
    void solve_poisson_u_p(VectorField& u_field, const VectorField& rhs, const VectorField& p_field);
    void solve_poisson_v_p(VectorField& v_field, const VectorField& rhs, const VectorField& p_field);

    // CG for pressure
    void solve_pressure_cg_custom();

    // CG helpers
    void calculation_r();          // r1 -= tau_cur * Ap  (in-place update)
    double scalar_product_rr();    // sets r_cur=(r1,r1), returns tau=(r1,r1)/(Ap,r1)

    // Output: saves field to output_dir/<filename>
    void save_field(const std::string& filename, const VectorField& data1, int boll);

    // Fields (public for easy inspection)
    VectorField u, v, p, u_tilde, v_tilde;
    VectorField fu, fv;   // RHS for u and v Poisson problems
    VectorField Ap;       // Result of Operator A
    VectorField r1;       // CG residual
    VectorField p_prev;   // Previous iterate
    VectorField F;        // Divergence: F = div(u_tilde, v_tilde)

    // F divergence auxiliary forces
    VectorField f11, f12, f22;

    // Bingham tensors
    VectorField tau11, tau12, tau22;
    VectorField gamma11, gamma12, gamma22;
    VectorField norma11, norma12;

    // Intermediate derivative and auxiliary fields
    VectorField du_dx, dv_dy; // Centers
    VectorField du_dy, dv_dx, D12; // Nodes
    VectorField t11, t12, t22; // Auxiliary tensors

    // FISTA/Uzawa auxiliary variables
    VectorField tau11_aux, tau12_aux, tau22_aux;
    VectorField tau11_aux_prev, tau12_aux_prev, tau22_aux_prev;
    double alpha_boost = 1.0, alpha_boost_prev = 1.0;

    // CG scalars
    double r_cur = 0.0, r_prev_val = 0.0;
    double tau_cur = 0.0, tau_prev = 0.0;
    double alpha_cur = 1.0, alpha_prev = 1.0;

    // Добавьте в секцию private или public
    void zero_field(VectorField& field);


private:
    LShapeConfig cfg;
    int nx_max, ny_max;

    // Node mappings
    std::map<std::pair<int,int>, int> u_node_to_k;
    std::vector<std::pair<int,int>> u_k_to_node;

    std::map<std::pair<int,int>, int> v_node_to_k;
    std::vector<std::pair<int,int>> v_k_to_node;

    // Pressure node mapping (cell centers)
    // Used by Operator_A to compute Bx*p and By*p
    // Pressure p lives at cell center (i-0.5, j-0.5): i=1..nx, j=1..ny

    // Eigen linear systems
    Eigen::SparseMatrix<double> Au, Av;
    Eigen::VectorXd b_u, b_v;
    Eigen::VectorXd sol_u, sol_v;
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver_u, solver_v;

    // Pre-allocated thread pool (2 workers: one for u, one for v)
    std::unique_ptr<ThreadPool> pool;

    // Domain helpers
    bool is_in_L_shape_u(int j, int i);
    bool is_in_L_shape_v(int j, int i);
    bool is_in_L_shape_p(int j, int i);
    std::string out_path(const std::string& name);
};
