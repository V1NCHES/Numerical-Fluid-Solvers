#include "PoissonSolver2D.h"
#include <iostream>
#include <iomanip>
#include <chrono>
/*
#define Ht  1e-1
#define H_2  Ht * Ht
#define DeltaT  0.1
#define SIGMA  H_2/(DeltaT*0.5)
#define diff_sigma2  SIGMA + 2.0
#define o_sum_sigma2 1/(SIGMA - 2.0)
*/
PoissonSolver2D::PoissonSolver2D(int nx, int ny, double lx, double ly)
{
    this->Nx = nx, this->Ny = ny, this->Lx = lx, this->Ly = ly;
    this->hx = Lx / Nx;
    this->hy = Ly / Ny;
    this->d_t = 0.01;
    this->h_y2 = hy * hy;
    double sigma = 2 * this->h_y2 / this->d_t;
    this->sm2 = sigma - 2;
    this->osp2 = 1/(sigma + 2);
}

// =============================================================================
// Вспомогательные функции для DST
// =============================================================================

double PoissonSolver2D::lambda_k_0(int k, double h, int N) {
    return 4.0 * std::sin(k * M_PI / (2.0 * N)) * std::sin(k * M_PI / (2.0 * N)) / (h * h);
}

double PoissonSolver2D::phi_k2(double** f_ij, int i, int k2, int N2) {
    double sum = 0.0;
    for (int j = 1; j < N2; j++) {
        sum += f_ij[i - 1][j - 1] * std::sin(k2 * M_PI * j / N2);
    }
    return sum;
}

double PoissonSolver2D::phi_k1_k2(double** f_ij, int k1, int N1, int k2, int N2) {
    double sum = 0.0;
    for (int i = 1; i < N1; i++) {
        sum += phi_k2(f_ij, i, k2, N2) * std::sin(k1 * M_PI * i / N1);
    }
    return sum;
}

double PoissonSolver2D::u_k2(double** f_ij, int i, double h1, int N1, int k2, double h2, int N2) {
    double sum = 0.0;
    for (int k1 = 1; k1 < N1; k1++) {
        double denominator = lambda_k_0(k1, h1, N1) + lambda_k_0(k2, h2, N2);
        if (denominator != 0.0) {
            sum += phi_k1_k2(f_ij, k1, N1, k2, N2) * std::sin(k1 * M_PI * i / N1) / denominator;
        }
    }
    return sum;
}

double PoissonSolver2D::u_ij(double** f_ij, int i, int j, double h1, double h2, int N1, int N2) {
    double sum = 0.0;
    for (int k2 = 1; k2 < N1; k2++) {
        sum += u_k2(f_ij, i, h1, N1, k2, h2, N2) * std::sin(k2 * M_PI * j / N2);
    }
    return sum;
}

// =============================================================================
// Метод дискретного синус-преобразования (DST)
// =============================================================================

void PoissonSolver2D::solve_DST(double** f_inner, double** u) 
{
    auto start_time = std::chrono::high_resolution_clock::now();

    for (int i = 1; i <= this->Nx - 1; i++)
    {
        for (int j = 1; j <= this->Ny - 1; j++) 
        {
            u[i - 1][j - 1] = u_ij(f_inner, i, j, hx, hy, Nx, Ny) * 4.0 / (Nx * Ny);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << " DST метод выполнен за " << duration.count() << " сек" << std::endl;

}


// =============================================================================
// Попеременно-треугольный метод
// =============================================================================
double PoissonSolver2D::U_ij(double** U, double** f, int i, int j)
{
    return (sm2 * U[i][j] + U[i + 1][j] + U[i - 1][j] + U[i][j + 1] + U[i][j - 1] + f[i][j] * h_y2) * osp2;
}

double PoissonSolver2D::Error(double** U, double** f, int i, int j)
{  
    return fabs(-4 * U[i][j] + U[i + 1][j] + U[i - 1][j] + U[i][j + 1] + U[i][j - 1] + f[i][j] * h_y2);
}
void PoissonSolver2D::solve_triangular_method(double** f, double** u)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    double ERROR = 0.0, ERROR1 = 0.0;
    double EPS = 1e-6;
    double Fmax = 0.0, Fmax1 = 0.0;
    Fmax = fabs(f[0][0] * h_y2);

    for (int j = 0; j < Ny - 1; ++j)
    {
        for (int i = 0; i < Nx - 1; ++i)
        {
            Fmax1 = fabs(f[i][j] * h_y2);
            if (Fmax1 > Fmax) Fmax = Fmax1;
        }
    }
    //std::cout << "Fmax = " << Fmax << "\n";

    int k = 0;
    do {

        for (int j = 1; j < Ny; ++j)
        {
            for (int i = 1; i < Nx; ++i)
            {
                u[i][j] = U_ij(u, f, i, j);
            }
        }

        for (int j = Ny - 1; j > 0; --j)
        {
            for (int i = Nx - 1; i > 0; --i)
            {
                u[i][j] = U_ij(u, f, i, j);
            }
        }

        ERROR = Error(u, f, 1, 1);

        for (int j = 1; j < Ny; ++j)
        {
            for (int i = 1; i < Nx; ++i)
            {
                ERROR1 = Error(u, f, i, j);
                if (ERROR1 > ERROR) ERROR = ERROR1;
            }
        }
        //*/
        //if (k % 1000 == 0) std::cout << " k = " << k << "\n";
        ++k;
    } while (ERROR > EPS * Fmax && k < 100000);

    //std::cout << "final k = " << k << "\n";
   
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "  Время: " << duration.count() << " сек" << std::endl;
}



// =============================================================================
// Вспомогательные функции
// =============================================================================

double PoissonSolver2D::f_test(double x, double y) {
    return -2.0 * (y * y - y + x * x - x);
}

std::vector<std::vector<double>> PoissonSolver2D::analytical_solution(int Nx, int Ny, double Lx, double Ly) {
    double hx = Lx / Nx;
    double hy = Ly / Ny;

    std::vector<std::vector<double>> analytical(Nx - 1, std::vector<double>(Ny - 1, 0.0));

    for (int i = 0; i < Nx - 1; i++) {
        for (int j = 0; j < Ny - 1; j++) {
            double x = (i + 1) * hx;
            double y = (j + 1) * hy;
            analytical[i][j] = (x * x - x) * (y * y - y);
        }
    }

    return analytical;
}

double PoissonSolver2D::compute_max_error(const std::vector<std::vector<double>>& numerical, const std::vector<std::vector<double>>& analytical)
{
    double max_error = 0.0;
    for (size_t i = 0; i < numerical.size(); i++)
    {
        for (size_t j = 0; j < numerical[i].size(); j++)
        {
            max_error = std::max(max_error, std::abs(numerical[i][j] - analytical[i][j]));
        }
    }
    return max_error;
}

double PoissonSolver2D::compute_mean_error(const std::vector<std::vector<double>>& numerical, const std::vector<std::vector<double>>& analytical)
{
    double sum_error = 0.0;
    int count = 0;

    for (size_t i = 0; i < numerical.size(); i++)
    {
        for (size_t j = 0; j < numerical[i].size(); j++)
        {
            sum_error += std::abs(numerical[i][j] - analytical[i][j]);
            count++;
        }
    }

    return (count > 0) ? sum_error / count : 0.0;
}

void PoissonSolver2D::write_to_file(const std::vector<std::vector<double>>& data, const std::string& filename, double Lx, double Ly)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
        return;
    }

    file << std::scientific << std::setprecision(6);

    size_t rows = data.size();
    size_t cols = (rows > 0) ? data[0].size() : 0;

    double dx = Lx / (rows + 1);
    double dy = Ly / (cols + 1);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            double x = (i + 1) * dx;
            double y = (j + 1) * dy;
            file << x << " " << y << " " << data[i][j] << std::endl;
        }
        file << std::endl;
    }
    /*
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                file << data[i][j] << " ";
            }
            file << std::endl;
        }
    */
    file.close();
    std::cout << "Данные записаны в файл: " << filename << std::endl;
}