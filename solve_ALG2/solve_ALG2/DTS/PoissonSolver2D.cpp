#include "PoissonSolver2D.h"
#include <iostream>
#include <iomanip>
#include <chrono>

PoissonSolver2D::PoissonSolver2D(int nx, int ny, double lx, double ly)
{
    this->Nx = nx, this->Ny = ny, this->Lx = lx, this->Ly = ly;
    this->hx = Lx / Nx;
    this->hy = Ly / Ny;
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
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "✓ DST метод выполнен за " << duration.count() << " мс" << std::endl;

    return u;
}

// =============================================================================
// Алгоритм Томаса
// =============================================================================

std::vector<double> PoissonSolver2D::thomas_algorithm_simple(const std::vector<double>& a, const std::vector<double>& b, const std::vector<double>& c, const std::vector<double>& d)
{
    int n = d.size();
    if (n == 0) return std::vector<double>();

    std::vector<double> a_copy = a;
    std::vector<double> b_copy = b;
    std::vector<double> c_copy = c;
    std::vector<double> d_copy = d;
    std::vector<double> x(n, 0.0);

    // Прямой ход
    for (int i = 1; i < n; i++) {
        double m = a_copy[i - 1] / b_copy[i - 1];
        b_copy[i] = b_copy[i] - m * c_copy[i - 1];
        d_copy[i] = d_copy[i] - m * d_copy[i - 1];
    }

    // Обратный ход
    x[n - 1] = d_copy[n - 1] / b_copy[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (d_copy[i] - c_copy[i] * x[i + 1]) / b_copy[i];
    }

    return x;
}

// =============================================================================
// Метод переменных направлений (ADI)
// =============================================================================

std::vector<std::vector<double>> PoissonSolver2D::solve_ADI(const std::vector<std::vector<double>>& f_inner, int max_iter, double tol) {
    auto start_time = std::chrono::high_resolution_clock::now();

    int Nx_inner = f_inner.size();
    int Ny_inner = f_inner[0].size();

    int total_points_x = Nx_inner + 2;
    int total_points_y = Ny_inner + 2;

    double hx_local = Lx / (total_points_x - 1);
    double hy_local = Ly / (total_points_y - 1);
    double hx2 = hx_local * hx_local;
    double hy2 = hy_local * hy_local;

    std::vector<std::vector<double>> u(Nx_inner, std::vector<double>(Ny_inner, 0.0));

    // Трехдиагональные матрицы
    std::vector<double> a_x(Nx_inner, -1.0 / hx2);
    std::vector<double> b_x(Nx_inner, 2.0 / hx2);
    std::vector<double> c_x(Nx_inner, -1.0 / hx2);

    std::vector<double> a_y(Ny_inner, -1.0 / hy2);
    std::vector<double> b_y(Ny_inner, 2.0 / hy2);
    std::vector<double> c_y(Ny_inner, -1.0 / hy2);

    for (int iteration = 0; iteration < max_iter; iteration++)
    {
        auto u_old = u;

        // Полушаг по X
        for (int j = 0; j < Ny_inner; j++) {
            std::vector<double> rhs(Nx_inner, 0.0);
            for (int i = 0; i < Nx_inner; i++) {
                double bottom = (j == 0) ? 0.0 : u[i][j - 1];
                double top = (j == Ny_inner - 1) ? 0.0 : u[i][j + 1];
                rhs[i] = -f_inner[i][j] - (bottom + top) / hy2;
            }

            std::vector<double> u_line = thomas_algorithm_simple(a_x, b_x, c_x, rhs);
            for (int i = 0; i < Nx_inner; i++) {
                u[i][j] = u_line[i];
            }
        }

        // Полушаг по Y
        for (int i = 0; i < Nx_inner; i++) {
            std::vector<double> rhs(Ny_inner, 0.0);
            for (int j = 0; j < Ny_inner; j++) {
                double left = (i == 0) ? 0.0 : u[i - 1][j];
                double right = (i == Nx_inner - 1) ? 0.0 : u[i + 1][j];
                rhs[j] = -f_inner[i][j] - (left + right) / hx2;
            }

            std::vector<double> u_column = thomas_algorithm_simple(a_y, b_y, c_y, rhs);
            for (int j = 0; j < Ny_inner; j++) {
                u[i][j] = u_column[j];
            }
        }

        // Проверка сходимости
        double residual = 0.0;
        for (int i = 0; i < Nx_inner; i++) {
            for (int j = 0; j < Ny_inner; j++) {
                residual = std::max(residual, std::abs(u[i][j] - u_old[i][j]));
            }
        }

        if (residual < tol) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            std::cout << "✓ ADI метод сошелся за " << iteration + 1 << " итераций" << std::endl;
            std::cout << "  Время: " << duration.count() << " мс, невязка: " << residual << std::endl;
            return u;
        }

        if (iteration % 100 == 0) {
            std::cout << "  Итерация " << iteration << ", невязка: " << residual << std::endl;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "✗ ADI метод не сошелся за " << max_iter << " итераций" << std::endl;
    std::cout << "  Время: " << duration.count() << " мс" << std::endl;

    return u;
}

// =============================================================================
// Попеременно-треугольный метод
// =============================================================================

std::vector<std::vector<double>> PoissonSolver2D::solve_triangular_method(const std::vector<std::vector<double>>& f_inner, int max_iter, double tol)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    int Nx_inner = f_inner.size();
    int Ny_inner = f_inner[0].size();

    int total_points_x = Nx_inner + 2;
    int total_points_y = Ny_inner + 2;

    double hx_local = Lx / (total_points_x - 1);
    double hy_local = Ly / (total_points_y - 1);
    double hx2 = hx_local * hx_local;
    double hy2 = hy_local * hy_local;

    double tau = 0.8 / (2.0 / hx2 + 2.0 / hy2);

    std::vector<std::vector<double>> u(Nx_inner, std::vector<double>(Ny_inner, 0.0));

    for (int iteration = 0; iteration < max_iter; iteration++)
    {
        auto u_old = u;

        // Прямой ход (черные узлы)
        for (int i = 0; i < Nx_inner; i++) {
            for (int j = 0; j < Ny_inner; j++) {
                if ((i + j) % 2 == 0) {
                    double left = (i == 0) ? 0.0 : u[i - 1][j];
                    double right = (i == Nx_inner - 1) ? 0.0 : u[i + 1][j];
                    double bottom = (j == 0) ? 0.0 : u[i][j - 1];
                    double top = (j == Ny_inner - 1) ? 0.0 : u[i][j + 1];

                    double laplacian_u = (left + right - 2 * u[i][j]) / hx2 + (bottom + top - 2 * u[i][j]) / hy2;
                    u[i][j] += tau * (f_inner[i][j] + laplacian_u);
                }
            }
        }

        // Обратный ход (красные узлы)
        for (int i = Nx_inner - 1; i >= 0; i--) {
            for (int j = Ny_inner - 1; j >= 0; j--) {
                if ((i + j) % 2 == 1) {
                    double left = (i == 0) ? 0.0 : u[i - 1][j];
                    double right = (i == Nx_inner - 1) ? 0.0 : u[i + 1][j];
                    double bottom = (j == 0) ? 0.0 : u[i][j - 1];
                    double top = (j == Ny_inner - 1) ? 0.0 : u[i][j + 1];

                    double laplacian_u = (left + right - 2 * u[i][j]) / hx2 + (bottom + top - 2 * u[i][j]) / hy2;
                    u[i][j] += tau * (f_inner[i][j] + laplacian_u);
                }
            }
        }

        // Проверка сходимости
        double residual = 0.0;
        for (int i = 0; i < Nx_inner; i++) {
            for (int j = 0; j < Ny_inner; j++) {
                residual = std::max(residual, std::abs(u[i][j] - u_old[i][j]));
            }
        }

        if (residual < tol) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            std::cout << "✓ Попеременно-треугольный метод сошелся за " << iteration + 1 << " итераций" << std::endl;
            std::cout << "  Время: " << duration.count() << " мс, невязка: " << residual << std::endl;
            return u;
        }

        if (iteration % 100 == 0) {
            std::cout << "  Итерация " << iteration << ", невязка: " << residual << std::endl;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "✗ Попеременно-треугольный метод не сошелся за " << max_iter << " итераций" << std::endl;
    std::cout << "  Время: " << duration.count() << " мс" << std::endl;

    return u;
}

// =============================================================================
// Шахматный метод (Red-Black)
// =============================================================================

std::vector<std::vector<double>> PoissonSolver2D::solve_checkerboard_method(const std::vector<std::vector<double>>& f_inner, int max_iter, double tol)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    int Nx_inner = f_inner.size();
    int Ny_inner = f_inner[0].size();

    int total_points_x = Nx_inner + 2;
    int total_points_y = Ny_inner + 2;

    double hx_local = Lx / (total_points_x - 1);
    double hy_local = Ly / (total_points_y - 1);
    double hx2 = hx_local * hx_local;
    double hy2 = hy_local * hy_local;

    std::vector<std::vector<double>> u(Nx_inner, std::vector<double>(Ny_inner, 0.0));

    for (int iteration = 0; iteration < max_iter; iteration++)
    {
        auto u_old = u;

        // Черные узлы
        for (int i = 0; i < Nx_inner; i++) {
            for (int j = 0; j < Ny_inner; j++) {
                if ((i + j) % 2 == 0) {
                    double left = (i == 0) ? 0.0 : u[i - 1][j];
                    double right = (i == Nx_inner - 1) ? 0.0 : u[i + 1][j];
                    double bottom = (j == 0) ? 0.0 : u[i][j - 1];
                    double top = (j == Ny_inner - 1) ? 0.0 : u[i][j + 1];

                    u[i][j] = (hx2 * hy2 * f_inner[i][j] +
                        hy2 * (left + right) +
                        hx2 * (bottom + top)) / (2.0 * (hx2 + hy2));
                }
            }
        }

        // Красные узлы
        for (int i = 0; i < Nx_inner; i++) {
            for (int j = 0; j < Ny_inner; j++) {
                if ((i + j) % 2 == 1) {
                    double left = (i == 0) ? 0.0 : u[i - 1][j];
                    double right = (i == Nx_inner - 1) ? 0.0 : u[i + 1][j];
                    double bottom = (j == 0) ? 0.0 : u[i][j - 1];
                    double top = (j == Ny_inner - 1) ? 0.0 : u[i][j + 1];

                    u[i][j] = (hx2 * hy2 * f_inner[i][j] +
                        hy2 * (left + right) +
                        hx2 * (bottom + top)) / (2.0 * (hx2 + hy2));
                }
            }
        }

        // Проверка сходимости
        double residual = 0.0;
        for (int i = 0; i < Nx_inner; i++) {
            for (int j = 0; j < Ny_inner; j++) {
                residual = std::max(residual, std::abs(u[i][j] - u_old[i][j]));
            }
        }

        if (residual < tol) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            std::cout << "✓ Шахматный метод сошелся за " << iteration + 1 << " итераций" << std::endl;
            std::cout << "  Время: " << duration.count() << " мс, невязка: " << residual << std::endl;
            return u;
        }

        if (iteration % 100 == 0) {
            std::cout << "  Итерация " << iteration << ", невязка: " << residual << std::endl;
        }
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "✗ Шахматный метод не сошелся за " << max_iter << " итераций" << std::endl;
    std::cout << "  Время: " << duration.count() << " мс" << std::endl;

    return u;
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

    file.close();
    std::cout << "✓ Данные записаны в файл: " << filename << std::endl;
}