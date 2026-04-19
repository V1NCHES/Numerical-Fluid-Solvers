#include <cmath>
#include <iostream>
#include "DST.h"
#include <chrono>



const double PI = 3.14159265358979323846;

double CorrectPoissonSolver2D::lamda_k_0(double k, double h, double N) {
    return 4.0 * pow(sin(k * PI / (2.0 * N)), 2) / (h * h);
}

double CorrectPoissonSolver2D::phi_k2(double** f_ij, int i, double k2, int N2) {
    double sum = 0.0;
    for (int j = 1; j < N2; j++) {
        if (i - 1 < 0 || j - 1 < 0) continue; // Защита от выхода за границы
        sum += f_ij[i - 1][j - 1] * sin(k2 * PI * j / N2);
    }
    return sum;
}

double CorrectPoissonSolver2D::phi_k1_k2(double** f_ij, double k1, int N1, double k2, int N2) {
    double sum = 0.0;
    for (int i = 1; i < N1; i++) {
        sum += phi_k2(f_ij, i, k2, N2) * sin(k1 * PI * i / N1);
    }
    return sum;
}

double CorrectPoissonSolver2D::u_k2(double** f_ij, int i, double h1, int N1, double k2, double h2, int N2) {
    double sum = 0.0;
    for (int k1 = 1; k1 < N1; k1++) {
        double numerator = phi_k1_k2(f_ij, k1, N1, k2, N2) * sin(k1 * PI * i / N1);
        double denominator = lamda_k_0(k1, h1, N1) + lamda_k_0(k2, h2, N2);
        if (denominator != 0.0) {
            sum += numerator / denominator;
        }
    }
    return sum;
}

double CorrectPoissonSolver2D::u_ij(double** f_ij, int i, int j, double h1, double h2, int N1, int N2) {
    double sum = 0.0;
    for (int k2 = 1; k2 < N1; k2++) {
        sum += u_k2(f_ij, i, h1, N1, k2, h2, N2) * sin(k2 * PI * j / N2);
    }
    return sum * 4.0 / (N1 * N2);
}

CorrectPoissonSolver2D::CorrectPoissonSolver2D(int N1, int N2, double L1, double L2) {
    this->N1 = N1;
    this->N2 = N2;
    this->L1 = L1;
    this->L2 = L2;
    this->h1 = L1 / N1;
    this->h2 = L2 / N2;

    std::cout << "Инициализация решателя: N1=" << N1 << ", N2=" << N2
        << ", h1=" << h1 << ", h2=" << h2 << std::endl;
}

void CorrectPoissonSolver2D::solve(double** f_inner, double** u) {
    std::cout << "Решение уравнения Пуассона..." << std::endl;
    std::cout << "Размер f_inner: " << N1 << "x" << N2 << std::endl;

    // Проверка входных данных
    if (f_inner == nullptr || u == nullptr) {
        std::cout << "ОШИБКА: nullptr в solve!" << std::endl;
        return;
    }

    for (int i = 1; i < N1; i++) {
        for (int j = 1; j < N2; j++) {
            if (i - 1 < N1 && j - 1 < N2) { // Проверка границ
                u[i - 1][j - 1] = u_ij(f_inner, i, j, h1, h2, N1, N2);
            }
        }
    }
    std::cout << "Решение уравнения Пуассона завершено" << std::endl;
}



