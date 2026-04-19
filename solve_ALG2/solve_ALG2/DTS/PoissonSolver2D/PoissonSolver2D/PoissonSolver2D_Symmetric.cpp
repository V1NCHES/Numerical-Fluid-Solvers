#include "PoissonSolver2D_Symmetric.h"
#include <iostream>
#include <iomanip>
#include <chrono>

// =============================================================================
// Конструктор
// =============================================================================
PoissonSolver2DSymmetric::PoissonSolver2DSymmetric(int nx, int ny, double lx, double ly, double d_t)
{
    this->Nx = nx; this->Ny = ny;
    this->Lx = lx;              // Длина по X
    this->Ly = ly;              // Длина по Y

    this->hx = this->Lx / (nx);
    this->hy = this->Ly / (ny);
    
    // Количество внутренних точек в половине области
    this->Nx_half = nx / 2 +1;
    this->Ny_half = ny / 2 + 1;

    this->d_t = d_t;
    this-> h_y2 = this->hy * this->hy;
    double sigma = 2.0 * this->h_y2 / d_t;
    this->sm2 = sigma - 2.0;
    double sp2 = sigma + 2.0;
    this->osp2 = 1 / sp2;

    // Проверка четности
    //if (nx % 2 != 0 || ny % 2 != 0) { std::cerr << "Ошибка: Количество внутренних точек должно быть четным для симметричного решения!" << std::endl;}

    // Шаг сетки: h = L / (N_inner + 1)
   

    std::cout << "Инициализация симметричного решателя:" << std::endl;
    std::cout << "  Полная область: " << Nx +1<< "x" << Ny +1 << "граница +  внутренних точек" << std::endl;
    std::cout << "  Половина области: " << Nx_half << "x" << Ny_half << " граница + внутренних точек" << std::endl;
    std::cout << "  Шаг сетки: hx = " << hx << ", hy = " << hy <<", d_t = " << d_t << std::endl;
}

PoissonSolver2DSymmetric::PoissonSolver2DSymmetric()
{
    this->Nx = 1; this->Ny = 1;  
    this->Lx = 1.0;              // Длина по X
    this->Ly = 1.0;              // Длина по Y
    this->hx = Lx / (Nx);
    this->hy = Ly / (Ny);


    this->d_t = 0.01;
    this->h_y2 = hy * hy;
    double sigma = 2.0 * h_y2 / d_t;
    this->sm2 = sigma - 2.0;
    double sp2 = sigma + 2.0;
    this->osp2 = 1 / sp2;
}

// =============================================================================
// Попеременно-треугольный метод для симеричной области
// =============================================================================
double PoissonSolver2DSymmetric::u_ij(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{
    return (this->sm2 * u[i][j] + u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] + f[i][j] * this->h_y2) * this->osp2;
}

double PoissonSolver2DSymmetric::u_ij_right(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{
    return ((this->sm2 +1)* u[i][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] + f[i][j] * this->h_y2) * this->osp2;
}

double PoissonSolver2DSymmetric::u_ij_above(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{
    return ((this->sm2 +1)* u[i][j] + u[i + 1][j] + u[i - 1][j] + u[i][j - 1] + f[i][j] * this->h_y2) * this->osp2;
}

double PoissonSolver2DSymmetric::u_ij_center(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{
    return ((this->sm2 +2)* u[i][j] + u[i - 1][j] + u[i][j - 1] + f[i][j] * this->h_y2) * this->osp2;
}

double PoissonSolver2DSymmetric::Error(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{
    return fabs(-4.0 * u[i][j] + u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] + f[i][j] * this->h_y2);
}

double PoissonSolver2DSymmetric::Error_right(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{
    return fabs(-3.0 * u[i][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] + f[i][j] * this->h_y2);    
}

double PoissonSolver2DSymmetric::Error_above(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{
    return fabs(-3.0 * u[i][j] + u[i + 1][j] + u[i - 1][j] + u[i][j - 1] + f[i][j] * this->h_y2);
}

double PoissonSolver2DSymmetric::Error_center(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f, int i, int j)
{
    return fabs(-2.0 * u[i][j] + u[i - 1][j] + u[i][j - 1] + f[i][j] * this->h_y2);
}
//////////////////////////////////////////////////////////////////////////////////

double PoissonSolver2DSymmetric::U_ij(double u_ij, double u_ip1j, double u_ijp1, double u_im1j, double u_ijm1, double f_ij)
{
    return (this->sm2 * u_ij + u_ip1j + u_im1j + u_ijp1 + u_ijm1 + f_ij * this->h_y2) * this->osp2;
}
double PoissonSolver2DSymmetric::U_ij_left(double u_ij, double u_ip1j, double u_ijp1,  double u_ijm1, double f_ij)
{
    return ((this->sm2 +1) * u_ij + u_ip1j + u_ijp1 + u_ijm1 + f_ij * this->h_y2) * this->osp2;
}
double PoissonSolver2DSymmetric::U_ij_down(double u_ij, double u_ip1j, double u_ijp1, double u_im1j, double f_ij)
{
    return ((this->sm2 +1)* u_ij + u_ip1j + u_im1j + u_ijp1  + f_ij * this->h_y2) * this->osp2;
}
double PoissonSolver2DSymmetric::U_ij_center(double u_ij, double u_ip1j, double u_ijp1,  double f_ij)
{
    return ((this->sm2 +2)* u_ij + u_ip1j + u_ijp1 + f_ij * this->h_y2) * this->osp2;
}


void PoissonSolver2DSymmetric::solve_triangular_method(std::vector<std::vector<double>>& u, std::vector<std::vector<double>>& f)
{
    //std::vector<std::vector<double>> Test(Nx + 1, std::vector<double>(Ny + 1, 0.0));
    double E = 0.0;
    double E1 = 0.0;
    double EPS = 1e-10;
    double Fmax = 0.0, Fmax1 = 0.0;
    Fmax = fabs(f[1][1] * h_y2);
    // int Nx_half, Ny_half;     // Количество точек в половине области
    for (int j = 1; j < Ny_half; ++j)
    {
        for (int i = 1; i < Nx_half; ++i)
        {
            Fmax1 = fabs(f[i][j] * h_y2);
            if (Fmax1 > Fmax) Fmax = Fmax1;
        }
    }
    std::cout << "Fmax = " << Fmax << "\n";

    
    int k = 0;
    do {

        for (int j = 1; j < Ny_half - 1; ++j)
        {
            for (int i = 1; i < Nx_half - 1; ++i)
            {
                u[i][j] = u_ij(u, f, i, j);
            }
            u[Nx_half - 1][j] = u_ij_right(u, f, Nx_half - 1, j);
        }

        for (int i = 1; i < Nx_half - 1; ++i)
        {
            u[i][Ny_half - 1] = u_ij_above(u, f,i, Ny_half - 1);
        }
        u[Nx_half - 1][Ny_half - 1] = u_ij_center(u, f, Nx_half - 1, Ny_half - 1);

        for (int j = 1; j < Ny_half - 1; ++j)
        {
            for (int i = 1; i < Nx_half - 1; ++i)
            {                  // ij      i+1j         ij+1         i-1j         ij-1
                u[i][j] = U_ij(u[i][j], u[i - 1][j], u[i][j - 1], u[i + 1][j], u[i][j + 1], f[i][j]);
            }
            u[Nx_half - 1][j] = U_ij_left(u[Nx_half - 1][j], u[Nx_half - 1 - 1][j], u[Nx_half - 1][j - 1], u[Nx_half - 1][j +1], f[Nx_half - 1][j]);
        }

        //std::cout << "\n\n";

        for (int i = 1; i < Nx_half - 1; ++i)
        {
            u[i][Ny_half-1] = U_ij_down(u[i][Ny_half-1], u[i - 1][Ny_half-1],  u[i + 1][Ny_half-1], u[i][Ny_half - 1 -1], f[i][Ny_half -1]);
        }

        u[Nx_half - 1][Ny_half - 1] = U_ij_center(u[Nx_half - 1][Ny_half - 1], u[Nx_half - 1 - 1][Ny_half - 1], u[Nx_half - 1][Ny_half - 1 - 1], f[Nx_half - 1][Ny_half - 1]);
  

        E = Error(u, f, 1, 1);
        for (int j = 1; j < Ny_half - 1; ++j)
        {
            for (int i = 1; i < Nx_half - 1; ++i)
            {
                E1 = Error(u, f, i, j);
                if (E1 > E) E = E1;
            }
        }

        for (int j = 1; j < Ny_half - 1; ++j)
        {
            E1 = Error_right(u, f, Nx_half - 1, j);
            if (E1 > E) E = E1;
        }

        for (int i = 1; i < Nx_half - 1; ++i)
        {
            E1 = Error_above(u, f, i, Ny_half - 1);   
            if (E1 > E) E = E1;
        }

        E1 = Error_center(u, f, Nx_half - 1, Ny_half - 1);
        if (E1 > E) E = E1;


            
        if (k % 100 == 0) std::cout << "E = "<< E<<" ; EPS*Fmax = "<< EPS * Fmax << " ; k = " << k << "\n";
        ++k;
    } while (E > EPS * Fmax && k < 10000); //1800

    std::cout << "final k = " << k << "\n";
    //return u;
}

// =============================================================================
// Попеременно-треугольный метод для области когда hx != hy
// =============================================================================



// =============================================================================
// Вспомогательные функции
// =============================================================================
void PoissonSolver2DSymmetric::symmetric_mapping(std::vector<std::vector<double>>& u)
{
    double cur = 0.0;
    std::cout << Nx_half << " ; " << Ny_half << "\n";


    for (int i = 1; i < Nx_half; ++i)
    {
        for (int j = 1; j < Ny_half; ++j)
        {
            cur = u[i][j];
            u[i][Ny - j] = cur;
            u[Nx-i][j] = cur;
            u[Nx-i ][Ny - j] = cur;
        }
    }

}
double PoissonSolver2DSymmetric::f_test(double x, double y) {
    return -2.0 * (y * y - y + x * x - x);
}

std::vector<std::vector<double>> PoissonSolver2DSymmetric::analytical_solution(int Nx, int Ny, double Lx, double Ly) {
    double hx = Lx / Nx;
    double hy = Ly / Ny;

    std::vector<std::vector<double>> analytical(Nx+1, std::vector<double>(Ny +1, 0.0));

    for (int i = 1; i < Nx; i++) {
        for (int j = 1; j < Ny; j++) {
            double x = (i) * hx;
            double y = (j) * hy;
            analytical[i][j] = (x * x - x) * (y * y - y);
        }
    }

    return analytical;
}

double PoissonSolver2DSymmetric::compute_max_error(const std::vector<std::vector<double>>& numerical, const std::vector<std::vector<double>>& analytical)
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

double PoissonSolver2DSymmetric::compute_mean_error(const std::vector<std::vector<double>>& numerical, const std::vector<std::vector<double>>& analytical)
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

void PoissonSolver2DSymmetric::write_to_file(const std::vector<std::vector<double>>& data, const std::string& filename)
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
            double x = (i)* this->hx;
            double y = (j) *this->hy;
            file << x << " " << y << " " << data[i][j] << std::endl;
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные записаны в файл: " << filename << std::endl;
}

void main_PoissonSolver2DSymmetric()
{
    // Установка русской локали
    setlocale(LC_ALL, "Russian");
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);

    std::cout << "Решение уравнения Пуассона для Симмеричного случая\n" << std::endl;
    std::cout << "==========================================================" << std::endl;

    // Параметры сетки
    int Nx = 101, Ny = 101;
    double L1 = 1.0, L2 = 1.0;
    double d_t = 0.001;
    // Создаем решатель
    PoissonSolver2DSymmetric solver(Nx, Ny, L1, L2, d_t);

    std::cout << "Параметры сетки:" << std::endl;
    std::cout << "  Nx = " << Nx << ", Ny = " << Ny << std::endl;
    std::cout << "  h1 = " << solver.get_hx() << ", h2 = " << solver.get_hy() << std::endl;
    std::cout << "  Общее число точек: " << solver.get_Nx() + 1 << " x " << solver.get_Ny() + 1 << std::endl;
    double x = 0.0;
    double y = 0.0;
    // Создаем правую часть
    std::vector<std::vector<double>> f(Nx + 1, std::vector<double>(Ny + 1, 0.0));
    std::vector<std::vector<double>> u(Nx + 1, std::vector<double>(Ny + 1, 0.0));
    for (int i = 1; i < Nx; ++i) {
        for (int j = 1; j < Ny; ++j) {
            x = (i)*solver.get_hx();
            y = (j)*solver.get_hy();
            f[i][j] = PoissonSolver2DSymmetric::f_test(x, y);
            //  std::cout << f[i][j] << " ";
        }
        // std::cout <<"\n";
    }


    // Аналитическое решение
    auto analytical = PoissonSolver2DSymmetric::analytical_solution(Nx, Ny, L1, L2);

    std::cout << "\nРешаем уравнение Пуассона различными методами:" << std::endl;


    // Решение попеременно-треугольным методом
    std::cout << "\n3. Попеременно-треугольный метод:" << std::endl;
    solver.solve_triangular_method(u, f);
    solver.symmetric_mapping(u);
    /*
    for (int i = 0; i < Nx + 1; ++i) {
        for (int j = 0; j < Ny + 1; ++j) {
            std::cout << f[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n\n^^^^^^^^^^^^^^^^^^^^^^^\n\n";
    for (int i = 1; i < Nx; ++i) {
        for (int j = 1; j < Ny; ++j) {
            std::cout << u[i][j] << " ";
        }
        std::cout <<"\n";
    }*/


    double max_error_triangular = PoissonSolver2DSymmetric::compute_max_error(u, analytical);
    double mean_error_triangular = PoissonSolver2DSymmetric::compute_mean_error(u, analytical);
    std::cout << "   Максимальная погрешность: " << max_error_triangular << std::endl;
    std::cout << "   Средняя погрешность: " << mean_error_triangular << std::endl;
    solver.write_to_file(u, "solution_triangular.txt");

    solver.write_to_file(analytical, "analytical_solution.txt");

    std::cout << "\n==========================================================" << std::endl;
    std::cout << "Все решения записаны в файлы:" << std::endl;
    std::cout << "==========================================================\n" << std::endl;


}