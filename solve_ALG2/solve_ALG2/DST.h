#pragma once

class CorrectPoissonSolver2D {
public:
    int N1;
    int N2;
    double L1;
    double L2;
    double h1;
    double h2;

    CorrectPoissonSolver2D(int N1, int N2, double L1, double L2);
    void solve(double** f_inner, double** u);

public:
    double lamda_k_0(double k, double h, double N);
    double phi_k2(double** f_ij, int i, double k2, int N2);
    double phi_k1_k2(double** f_ij, double k1, int N1, double k2, int N2);
    double u_k2(double** f_ij, int i, double h1, int N1, double k2, double h2, int N2);
    double u_ij(double** f_ij, int i, int j, double h1, double h2, int N1, int N2);
};