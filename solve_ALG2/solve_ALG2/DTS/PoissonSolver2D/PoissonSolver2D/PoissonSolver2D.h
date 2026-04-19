#ifndef POISSONSOLVER2D_H
#define POISSONSOLVER2D_H

#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <locale>
#include <windows.h>
#define M_PI acos(-1)
//#define max_iter 10000
//#define tol 1e-6


class PoissonSolver2D {
private:
    int Nx, Ny;           // Количество интервалов
    double Lx, Ly;        // Размеры области
    double hx, hy;        // Шаги сетки
    double d_t;
    // Вспомогательные функции
    double lambda_k_0(int k, double h, int N);
    double phi_k2(const std::vector<std::vector<double>>& f_ij, int i, int k2, int N2);
    double phi_k1_k2(const std::vector<std::vector<double>>& f_ij, int k1, int N1, int k2, int N2);
    double u_k2(const std::vector<std::vector<double>>& f_ij, int i, double h1, int N1, int k2, double h2, int N2);
    double u_ij(const std::vector<std::vector<double>>& f_ij, int i, int j, double h1, double h2, int N1, int N2);

 

public:
    PoissonSolver2D(int nx, int ny, double lx, double ly);
    PoissonSolver2D(double Hx, double Hy, double lx, double ly);

    // Основные методы решения
    std::vector<std::vector<double>> solve_DST(const std::vector<std::vector<double>>& f_inner);
    std::vector<std::vector<double>> solve_triangular_method(std::vector<std::vector<double>>& u,  std::vector<std::vector<double>>& f_inner);
    double U_ij(std::vector<std::vector<double>>& U, std::vector<std::vector<double>>& f, int i, int j);
    double Error(std::vector<std::vector<double>>& U, std::vector<std::vector<double>>& f, int i, int j);
    // Вспомогательные функции
    static double f_test(double x, double y);
    static std::vector<std::vector<double>> analytical_solution(int Nx, int Ny, double Lx, double Ly);

    // Вычисление ошибки
    static double compute_max_error( std::vector<std::vector<double>>& numerical,  std::vector<std::vector<double>>& analytical);
    static double compute_mean_error( std::vector<std::vector<double>>& numerical,  std::vector<std::vector<double>>& analytical);

    // Запись результатов в файл
    static void write_to_file(const std::vector<std::vector<double>>& data, const std::string& filename, double Lx, double Ly);

    // Геттеры
    int get_Nx() const { return Nx; }
    int get_Ny() const { return Ny; }
    double get_hx() const { return hx; }
    double get_hy() const { return hy; }

    
};
void main_PoissonSolver2D();

#endif

