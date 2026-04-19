#ifndef POISSONSOLVER2D_SYMMETRIC_H
#define POISSONSOLVER2D_SYMMETRIC_H

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

class PoissonSolver2DSymmetric {
private:
    int Nx, Ny; // Количесво квадратиков
    int Nx_half, Ny_half;     // Количество точек в половине области
    double Lx, Ly;            // Размеры области
    double hx, hy;            // Шаги сетки

    // парамеры Для ПТМ
    double d_t;                // Шаг по времени
    double h_y2;                // Шаг^2 по y; sigma = 2*h_y2/d_t ;
    double sm2;                  // sigma + 2                           
    double osp2;                // 1/(sigma-2)

public:
    // Конструктор принимает количество Количесво квадратиков
    PoissonSolver2DSymmetric(int nx, int ny, double lx, double ly, double d_t);
    PoissonSolver2DSymmetric();

    double u_ij(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j);
    double u_ij_right(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j);
    double u_ij_above(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j);
    double u_ij_center(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j);

    double Error(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j);
    double Error_right(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j);
    double Error_above(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j);
    double Error_center(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j);


    double U_ij(double u_ij, double u_ip1j, double u_ijp1, double u_im1j, double u_ijm1, double f_ij);
    double U_ij_left(double u_ij, double u_ip1j, double u_ijp1, double u_ijm1, double f_ij);
    double U_ij_down(double u_ij, double u_ip1j, double u_ijp1, double u_im1j, double f_ij);
    double U_ij_center(double u_ij, double u_ip1j, double u_ijp1, double f_ij);



    void solve_triangular_method(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f);


    void symmetric_mapping(std::vector<std::vector<double>>& u);
    // Вспомогательные функции
    static double f_test(double x, double y);
    static std::vector<std::vector<double>> analytical_solution(int Nx, int Ny, double Lx, double Ly);

    // Вычисление ошибки
    static double compute_max_error(const std::vector<std::vector<double>>& numerical, const std::vector<std::vector<double>>& analytical);
   static double compute_mean_error(const std::vector<std::vector<double>>& numerical, const std::vector<std::vector<double>>& analytical);

    // Запись результатов в файл
    void write_to_file(const std::vector<std::vector<double>>& data, const std::string& filename);


    // Геттеры
    int get_Nx_half() const { return Nx_half; }
    int get_Ny_half() const { return Ny_half; }
    double get_hx() const { return hx; }
    double get_hy() const { return hy; }
    double get_Nx() const { return Nx; }
    double get_Ny() const { return Ny; }

    
};
void main_PoissonSolver2DSymmetric();
#endif