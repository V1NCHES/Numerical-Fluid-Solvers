#include "PoissonSolver2D.h"
#include <iostream>
#include <iomanip>
#include <chrono>

//const double h_t = 0.001;

//const double h_t2 = h_t* h_t;
//const double sigma = 2 * h_t2 / d_t;
//const double sm2 = sigma - 2;
//const double sp2 = sigma + 2;
//const double osp2 = 1/ sp2;

PoissonSolver2D::PoissonSolver2D(int nx, int ny, double lx, double ly)
{
    this->Nx = nx, this->Ny = ny, this->Lx = lx, this->Ly = ly;
    this->hx = Lx / Nx;
    this->hy = Ly / Ny;
    this-> d_t = 0.001;
    std::cout << "hx = " << hx << "\n";
}

PoissonSolver2D::PoissonSolver2D(double Hx, double Hy, double lx, double ly)
{
    this->Lx = lx, this->Ly = ly;
    this->hx = Hx;
    this->hy = Hy;
    this->d_t = 0.001;
    
    this->Nx = Lx/hx, this->Ny = Ly/hy;
    std::cout << "Nx = " << Nx << "; Ny = "<< Ny << "\n";
}



// =============================================================================
// Вспомогательные функции для DST
// =============================================================================

double PoissonSolver2D::lambda_k_0(int k, double h, int N) {
    return 4.0 * std::sin(k * M_PI / (2.0 * N)) * std::sin(k * M_PI / (2.0 * N)) / (h * h);
}

double PoissonSolver2D::phi_k2(const std::vector<std::vector<double>>& f_ij, int i, int k2, int N2) {
    double sum = 0.0;
    for (int j = 1; j < N2; j++) {
        sum += f_ij[i - 1][j - 1] * std::sin(k2 * M_PI * j / N2);
    }
    return sum;
}

double PoissonSolver2D::phi_k1_k2(const std::vector<std::vector<double>>& f_ij, int k1, int N1, int k2, int N2) {
    double sum = 0.0;
    for (int i = 1; i < N1; i++) {
        sum += phi_k2(f_ij, i, k2, N2) * std::sin(k1 * M_PI * i / N1);
    }
    return sum;
}

double PoissonSolver2D::u_k2(const std::vector<std::vector<double>>& f_ij, int i, double h1, int N1, int k2, double h2, int N2) {
    double sum = 0.0;
    for (int k1 = 1; k1 < N1; k1++) {
        double denominator = lambda_k_0(k1, h1, N1) + lambda_k_0(k2, h2, N2);
        if (denominator != 0.0) {
            sum += phi_k1_k2(f_ij, k1, N1, k2, N2) * std::sin(k1 * M_PI * i / N1) / denominator;
        }
    }
    return sum;
}

double PoissonSolver2D::u_ij(const std::vector<std::vector<double>>& f_ij, int i, int j, double h1, double h2, int N1, int N2) {
    double sum = 0.0;
    for (int k2 = 1; k2 < N1; k2++) {
        sum += u_k2(f_ij, i, h1, N1, k2, h2, N2) * std::sin(k2 * M_PI * j / N2);
    }
    return sum;
}

// =============================================================================
// Метод дискретного синус-преобразования (DST)
// =============================================================================

std::vector<std::vector<double>> PoissonSolver2D::solve_DST(const std::vector<std::vector<double>>& f_inner) {
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


// =============================================================================
// Попеременно-треугольный метод
// =============================================================================
double PoissonSolver2D::U_ij(std::vector<std::vector<double>>& U, std::vector<std::vector<double>>& f, int i, int j)
{
    double h_y2 = hy*hy;
    const double sigma = 2 * h_y2 / d_t;
    const double sm2 = sigma - 2;
    const double sp2 = sigma + 2;
    const double osp2 = 1 / sp2;
    return (sm2 * U[i][j] + U[i+1][j] + U[i-1][j] + U[i][j+1] + U[i][j-1] + f[i][j]* h_y2) * osp2;
}

double PoissonSolver2D::Error(std::vector<std::vector<double>>& U, std::vector<std::vector<double>>& f, int i, int j)
{
    double h_y2 = hy * hy;
    return fabs(-4 *U[i][j] +  U[i + 1][j] + U[i - 1][j] + U[i][j + 1] + U[i][j - 1] + f[i][j] * h_y2);
}

std::vector<std::vector<double>> PoissonSolver2D::solve_triangular_method(std::vector<std::vector<double>>&u,  std::vector<std::vector<double>>& f_inner)
{
    std::vector<std::vector<double>> f(Nx + 1, std::vector<double>(Ny + 1, 0.0));
    std::vector<std::vector<double>> U(Nx + 1, std::vector<double>(Ny + 1, 0.0));
    double h_y2 = hy * hy;
    double Error0 = 0.0, ERROR1 = 0.0;
    double EPS = 1e-10;
    double Fmax = 0.0, Fmax1 = 0.0;
    Fmax = fabs(f_inner[1][1] * h_y2);
    
    for (int j = 1; j < Ny-1 ; ++j)
    {
        for (int i = 1; i < Nx-1; ++i)
        {
            Fmax1 = fabs(f_inner[i][j] * h_y2);
            if (Fmax1 > Fmax) Fmax = Fmax1;
        }
    }
    std::cout <<"Fmax = "<< Fmax<<"\n";

   
    for (int i = 0; i < Nx-1; ++i)
    {
        for (int j = 0; j < Ny-1; ++j)
        {
            U[i+1][j+1] = u[i][j];
            f[i+1][j+1] = f_inner[i][j];
        }
    }
    /*
    for (int i = 0; i < Nx+1 ; ++i)
    {
        for (int j = 0; j < Ny+1 ; ++j)
        {
            std::cout << U[i][j] << " ";
        }
        std::cout << "\n";
    }
    */
    int k = 0;
    do {
       
        for (int j = 1; j < Ny; ++j)
        {
            for (int i = 1; i < Nx; ++i)
            {
                U[i][j] = U_ij(U, f, i, j);              
            }
        }
       
        for (int j = Ny-1; j > 0; --j)
        {
            for (int i = Nx-1; i > 0; --i)
            {
                U[i][j] = U_ij(U, f, i, j);
            }
        }
        /*
        for (int j = Ny-1; j > 0; --j)
        {
            for (int i = 1; i < Nx; ++i)
            {
                U[i][j] = U_ij(U, f, i, j);
            }
        }

        for (int j = 0; j < Ny; ++j)
        {
            for (int i = Nx-1; i > 0; --i)
            {
                U[i][j] = U_ij(U, f, i, j);
            }
        }*/
        
        Error0 = Error(U, f, 1, 1);
        
        for (int j = 1; j < Ny; ++j)
        {
            for (int i = 1; i < Nx; ++i)
            {
                ERROR1 = Error(U, f, i, j);
                if (ERROR1 > Error0) Error0 = ERROR1;
            }
        }
        //*/
        if (k % 1000 == 0) std::cout << " k = " << k << "\n";
        ++k;
    } while (Error0 > EPS * Fmax && k < 100000);
   /*
    for (int i = 0; i < Nx + 1; ++i)
    {
        for (int j = 0; j < Ny + 1; ++j)
        {
            std::cout << U[i][j] << " ";
        }
        std::cout << "\n";
    }*/

    for (int i = 0; i < Nx-1; ++i)
    {
        for (int j = 0; j < Ny-1; ++j)
        {
            u[i][j] = U[i + 1][j + 1];
            //std::cout << u[i][j] << " ";
        }
       // std::cout << "\n";
    }
    //*/
    std::cout << "final k = " << k << "\n";
    return u;
}



// =============================================================================
// Вспомогательные функции
// =============================================================================

double PoissonSolver2D::f_test(double x, double y) {
    //return -2.0 * (y * y - y + x * x - x);
    return 2.0 * M_PI * M_PI * sin(M_PI * x) * sin(M_PI * x);
}

std::vector<std::vector<double>> PoissonSolver2D::analytical_solution(int Nx, int Ny, double Lx, double Ly) {
    double hx = Lx / Nx;
    double hy = Ly / Ny;

    std::vector<std::vector<double>> analytical(Nx - 1, std::vector<double>(Ny - 1, 0.0));

    for (int i = 0; i < Nx - 1; i++) {
        for (int j = 0; j < Ny - 1; j++) {
            double x = (i + 1) * hx;
            double y = (j + 1) * hy;
            //analytical[i][j] = (x * x - x) * (y * y - y);
            analytical[i][j] = sin(M_PI * x) * sin(M_PI * y);
        }
    }

    return analytical;
}

double PoissonSolver2D::compute_max_error( std::vector<std::vector<double>>& numerical,  std::vector<std::vector<double>>& analytical)
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

double PoissonSolver2D::compute_mean_error( std::vector<std::vector<double>>& numerical,  std::vector<std::vector<double>>& analytical)
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

    int rows = data.size();
    int cols = (rows > 0) ? data[0].size() : 0;

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



void main_PoissonSolver2D()
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
   // int N1 = 100, N2 = 200;
    double L1 = 1.0, L2 = 1.0;
    double Hx = 0.005, Hy = 0.005;

    // Создаем решатель
    //PoissonSolver2D solver(N1, N2, L1, L2);
    PoissonSolver2D solver(Hx, Hy, L1, L2);
    int N1 = solver.get_Nx(), N2 = solver.get_Ny();
    std::cout << "Параметры сетки:" << std::endl;
    std::cout << "  N1 = " << N1 << ", N2 = " << N2 << std::endl;
    std::cout << "  h1 = " << solver.get_hx() << ", h2 = " << solver.get_hy() << std::endl;
    std::cout << "  Общее число точек: " << solver.get_Nx() + 1 << " x " << solver.get_Ny() + 1 << std::endl;

    // Создаем правую часть
    std::vector<std::vector<double>> f(N1-1, std::vector<double>(N2-1, 0.0));
    std::vector<std::vector<double>> u(N1-1, std::vector<double>(N2-1, 0.0));
    for (int i = 0; i < N1-1; i++) {
        for (int j = 0; j < N2-1; j++) {
            double x = (i + 1) * solver.get_hx();
            double y = (j + 1) * solver.get_hy();
            f[i][j] = PoissonSolver2D::f_test(x, y);
        }
    }
    // Аналитическое решение
    auto analytical = PoissonSolver2D::analytical_solution(N1, N2, L1, L2);
    std::cout << "\nРешаем уравнение Пуассона различными методами:" << std::endl;
    // Решение DST методом
    /*
    std::cout << "\n1. Метод дискретного синус-преобразования (DST):" << std::endl;
    auto solution_DST = solver.solve_DST(f);
    double max_error_DST = PoissonSolver2D::compute_max_error(solution_DST, analytical);
    double mean_error_DST = PoissonSolver2D::compute_mean_error(solution_DST, analytical);
    std::cout << "   Максимальная погрешность: " << max_error_DST << std::endl;
    std::cout << "   Средняя погрешность: " << mean_error_DST << std::endl;
    PoissonSolver2D::write_to_file(solution_DST, "solution_DST.txt", L1, L2);
    */
    // Решение попеременно-треугольным методом
   
    std::cout << "\n3. Попеременно-треугольный метод:" << std::endl;
    auto solution_triangular = solver.solve_triangular_method(u, f);
    double max_error_triangular = PoissonSolver2D::compute_max_error(solution_triangular, analytical);
    double mean_error_triangular = PoissonSolver2D::compute_mean_error(solution_triangular, analytical);
    std::cout << "   Максимальная погрешность: " << max_error_triangular << std::endl;
    std::cout << "   Средняя погрешность: " << mean_error_triangular << std::endl;
    PoissonSolver2D::write_to_file(solution_triangular, "solution_triangular.txt", L1, L2);

    PoissonSolver2D::write_to_file(analytical, "analytical_solution.txt", L1, L2);
    //*/
    std::cout << "\n==========================================================" << std::endl;
    std::cout << "Все решения записаны в файлы:" << std::endl;
    std::cout << "==========================================================\n" << std::endl;
}
