#include "PoissonSolver2D_hx_hy.h"
#include <iostream>
#include <iomanip>
#include <chrono>

//const double h_t = 0.001;

//const double h_t2 = h_t* h_t;
//const double sigma = 2 * h_t2 / d_t;
//const double sm2 = sigma - 2;
//const double sp2 = sigma + 2;
//const double osp2 = 1/ sp2;

PoissonSolver2D_hx_hy::PoissonSolver2D_hx_hy(int nx, int ny, double lx, double ly, double d_t)
{
    this->Nx = nx; this->Ny = ny;
    this->Lx = lx;              // Длина по X
    this->Ly = ly;              // Длина по Y

    this->hx = this->Lx / (nx);
    this->hy = this->Ly / (ny);

    this->hx2 = hx * hx, this->hy2 = hy * hy;

    this->d_t = d_t;
    this->ax = 2 * hx2 / d_t;
    this->ay = 2 * hy2 / d_t;

    this->r = hx2 / hy2;
    this->q = 1 / r;

    u_prev = std::vector<std::vector<double>>(Nx + 1, std::vector<double>(Ny + 1, 0.0));
   
    
    this->o_ax2 = 1/(ax+2), this->o_ay2 = 1 / (ay + 2);
    std::cout << "Инициализация симметричного решателя:" << std::endl;
    std::cout << "  Полная область: " << Nx + 1 << "x" << Ny + 1 << "граница +  внутренних точек" << std::endl;
    std::cout << "  Шаг сетки: hx = " << hx << ", hy = " << hy << ", d_t = " << d_t << std::endl;
    std::cout << "  Парамеры: o_ax2 = " << o_ax2 << ", o_ay2 = " << o_ay2 << ", r = " << r << ", q = " << q<<std::endl;
}

PoissonSolver2D_hx_hy::PoissonSolver2D_hx_hy()
{
    this->Nx = 1; this->Ny = 1;
    this->Lx = 1.0;              // Длина по X
    this->Ly = 1.0;              // Длина по Y
    this->hx = Lx / (Nx);
    this->hy = Ly / (Ny);

    this->d_t = 0.01;
    this->ax = 2*hx*hx/d_t;
    this->ay = 2 * hy * hy / d_t;

    this->r = hx*hx/(hy*hy);
    this->q = 1/r;
    u_prev = std::vector<std::vector<double>>(Nx + 1, std::vector<double>(Ny + 1, 0.0));
}

// =============================================================================
// Вспомогательные функции для DST
// =============================================================================
/*
double PoissonSolver2D_hx_hy::lambda_k_0(int k, double h, int N) {
    return 4.0 * std::sin(k * M_PI / (2.0 * N)) * std::sin(k * M_PI / (2.0 * N)) / (h * h);
}

double PoissonSolver2D_hx_hy::phi_k2(const std::vector<std::vector<double>>& f_ij, int i, int k2, int Ny) {
    double sum = 0.0;
    for (int j = 1; j < Ny; j++) {
        sum += f_ij[i - 1][j - 1] * std::sin(k2 * M_PI * j / Ny);
    }
    return sum;
}

double PoissonSolver2D_hx_hy::phi_k1_k2(const std::vector<std::vector<double>>& f_ij, int k1, int Nx, int k2, int Ny) {
    double sum = 0.0;
    for (int i = 1; i < Nx; i++) {
        sum += phi_k2(f_ij, i, k2, Ny) * std::sin(k1 * M_PI * i / Nx);
    }
    return sum;
}

double PoissonSolver2D_hx_hy::u_k2(const std::vector<std::vector<double>>& f_ij, int i, double h1, int Nx, int k2, double h2, int Ny) {
    double sum = 0.0;
    for (int k1 = 1; k1 < Nx; k1++) {
        double denominator = lambda_k_0(k1, h1, Nx) + lambda_k_0(k2, h2, Ny);
        if (denominator != 0.0) {
            sum += phi_k1_k2(f_ij, k1, Nx, k2, Ny) * std::sin(k1 * M_PI * i / Nx) / denominator;
        }
    }
    return sum;
}

double PoissonSolver2D_hx_hy::u_ij(const std::vector<std::vector<double>>& f_ij, int i, int j, double h1, double h2, int Nx, int Ny) {
    double sum = 0.0;
    for (int k2 = 1; k2 < Nx; k2++) {
        sum += u_k2(f_ij, i, h1, Nx, k2, h2, Ny) * std::sin(k2 * M_PI * j / Ny);
    }
    return sum;
}

// =============================================================================
// Метод дискретного синус-преобразования (DST)
// =============================================================================

std::vector<std::vector<double>> PoissonSolver2D_hx_hy::solve_DST(const std::vector<std::vector<double>>& f_inner) {
    auto start_time = std::chrono::high_resolution_clock::now();

    int inner_Nx = f_inner.size();
    int inner_Ny = f_inner[0].size();
    std::vector<std::vector<double>> u(inner_Nx, std::vector<double>(inner_Ny, 0.0));

    for (int i = 1; i <= inner_Nx; i++) {
        for (int j = 1; j <= inner_Ny; j++) {
            u[i - 1][j - 1] = u_ij(f_inner, i, j, hx, hy, Nx, Ny) * 4.0 / (Nx * Ny);
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << "DST метод выполнен за " << duration.count() << " сек" << std::endl;

    return u;
}
*/

// =============================================================================
// Попеременно-треугольный метод
// =============================================================================
double PoissonSolver2D_hx_hy::U_ij_x(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{
    
    return (u[i - 1][j] + u[i][j - 1] + ax* u_prev[i][j]  + r*(u_prev[i][j+1] + u_prev[i][j - 1] - 2 * u_prev[i][j])  + f[i][j] * hx2) * o_ax2;
}

double PoissonSolver2D_hx_hy::U_ij_y(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{

    return (u[i + 1][j] + u[i][j + 1] + ay * u_prev[i][j] + q * (u_prev[i+1][j] + u_prev[i-1][j] -  2 * u_prev[i][j]) + f[i][j] * hy2) * o_ay2;
}

double PoissonSolver2D_hx_hy::Error(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{ 
    return fabs(u[i + 1][j] + u[i - 1][j] - 2* u[i][j] + r*(u[i][j + 1] + u[i][j - 1] - 2 * u[i][j]) + f[i][j] * hx2);
}

void PoissonSolver2D_hx_hy::solve_triangular_method(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f)
{
   
    double Error0 = 0.0, ERROR1 = 0.0;
    double EPS = 1e-8;
    double Fmax = 0.0, Fmax1 = 0.0;
    Fmax = fabs(f[1][1] * hx*hx);

    for (int j = 1; j < Ny + 1; ++j)
    {
        for (int i = 1; i < Nx + 1; ++i)
        {
            Fmax1 = fabs(f[i][j] * hx * hx);
            if (Fmax1 > Fmax) Fmax = Fmax1;
        }
    }
    std::cout << "Fmax = " << Fmax << "\n";

    int k = 0;
    do {
        u_prev = u;
        for (int j = 1; j < Ny; ++j)
        {
            for (int i = 1; i < Nx; ++i)
            {
                u[i][j] = U_ij_x(u, f, i, j);
            }
        }
        u_prev = u;
        for (int j = Ny-1; j > 0; --j)
        {
            for (int i = Nx-1; i > 0; --i)
            {
                u[i][j] = U_ij_y(u, f, i, j);
            }
        }
       

        Error0 = Error(u, f, 1, 1);

        for (int j = 1; j < Ny; ++j)
        {
            for (int i = 1; i < Nx; ++i)
            {
                ERROR1 = Error(u, f, i, j);
                if (ERROR1 > Error0) Error0 = ERROR1;
            }
        }
        //*/
        if (k % 1000 == 0) std::cout << " k = " << k << "\n";
        ++k;
    } while (Error0 > EPS * Fmax && k < 100000);
    
    std::cout << "final k = " << k << "\n";
}



// =============================================================================
// Вспомогательные функции
// =============================================================================

double PoissonSolver2D_hx_hy::f_test(double x, double y) {
    return -2.0 * (y * y - y + x * x - x);
}

std::vector<std::vector<double>> PoissonSolver2D_hx_hy::analytical_solution(int Nx, int Ny, double Lx, double Ly) {
    double hx = Lx / Nx;
    double hy = Ly / Ny;

    std::vector<std::vector<double>> analytical(Nx + 1, std::vector<double>(Ny + 1, 0.0));

    for (int i = 1; i < Nx; i++) {
        for (int j = 1; j < Ny; j++) {
            double x = (i)*hx;
            double y = (j)*hy;
            analytical[i][j] = (x * x - x) * (y * y - y);
        }
    }

    return analytical;
}

double PoissonSolver2D_hx_hy::compute_max_error(const std::vector<std::vector<double>>& numerical, const std::vector<std::vector<double>>& analytical)
{
    double max_error = 0.0;
    for (size_t i = 0; i < numerical.size(); i++)
    {
        for (size_t j = 0; j < numerical[i].size(); j++)
        {
            max_error = max(max_error, std::abs(numerical[i][j] - analytical[i][j]));
        }
    }
    return max_error;
}

double PoissonSolver2D_hx_hy::compute_mean_error(const std::vector<std::vector<double>>& numerical, const std::vector<std::vector<double>>& analytical)
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

void PoissonSolver2D_hx_hy::write_to_file(const std::vector<std::vector<double>>& data, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
        return;
    }

    file << std::scientific << std::setprecision(6);

    for (int i = 0; i < Nx + 1; i++)
    {
        for (int j = 0; j < Ny + 1; j++)
        {
            double x = (i)*this->hx;
            double y = (j)*this->hy;
            file << x << " " << y << " " << data[i][j] << std::endl;
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные записаны в файл: " << filename << std::endl;
}



void main_PoissonSolver2D_hx_hy()
{

    // Установка русской локали
    setlocale(LC_ALL, "Russian");
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);

    std::cout << "Решение уравнения Пуассона\n" << std::endl;
    std::cout << "==========================================================" << std::endl;

    // Параметры сетки
    int Nx = 10, Ny = 10;
    double L1 = 1.0, L2 = 1.0;
    double d_t = 0.01;
    // Создаем решатель
    PoissonSolver2D_hx_hy solver(Nx, Ny, L1, L2, d_t);

    std::cout << "Параметры сетки:" << std::endl;
    std::cout << "  Nx = " << Nx << ", Ny = " << Ny << std::endl;
    std::cout << "  h1 = " << solver.get_hx() << ", h2 = " << solver.get_hy() << std::endl;
    std::cout << "  Общее число точек: " << solver.get_Nx() + 1 << " x " << solver.get_Ny() + 1 << std::endl;

    // Создаем правую часть
    std::vector<std::vector<double>> f(Nx + 1, std::vector<double>(Ny + 1, 0.0));
    std::vector<std::vector<double>> u(Nx + 1, std::vector<double>(Ny + 1, 0.0));
    for (int i = 1; i < Nx ; i++) {
        for (int j = 1; j < Ny; j++) {
            double x = (i) * solver.get_hx();
            double y = (j) * solver.get_hy();
            f[i][j] = PoissonSolver2D_hx_hy::f_test(x, y);
        }
    }
    // Аналитическое решение
    auto analytical = PoissonSolver2D_hx_hy::analytical_solution(Nx, Ny, L1, L2);
    std::cout << "\nРешаем уравнение Пуассона различными методами:" << std::endl;

    // Решение попеременно-треугольным методом
    std::cout << "\n3. Попеременно-треугольный метод:" << std::endl;
    solver.solve_triangular_method(u, f);

   /*
   for (int i = 0; i < Nx + 1; ++i) {
       for (int j = 0; j < Ny + 1; ++j) {
           std::cout << f[i][j] << " ";
       }
       std::cout << "\n";
   }*/
   std::cout << "\n\n^^^^^^^^^^^^^^^^^^^^^^^\n\n";
   for (int i = 0; i < Nx+1; ++i) {
       for (int j = 0; j < Ny+1; ++j) {
           std::cout << u[i][j] << " ";
       }
       std::cout <<"\n";
   }//*/


    double max_error_triangular = PoissonSolver2D_hx_hy::compute_max_error(u, analytical);
    double mean_error_triangular = PoissonSolver2D_hx_hy::compute_mean_error(u, analytical);
    std::cout << "   Максимальная погрешность: " << max_error_triangular << std::endl;
    std::cout << "   Средняя погрешность: " << mean_error_triangular << std::endl;
    solver.write_to_file(u, "solution_triangular.txt");

    solver.write_to_file(analytical, "analytical_solution.txt");

    std::cout << "\n==========================================================" << std::endl;
    std::cout << "Все решения записаны в файлы:" << std::endl;
    std::cout << "==========================================================\n" << std::endl;
}
