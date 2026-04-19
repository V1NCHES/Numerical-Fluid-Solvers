#include "PoissonSolver2D.h"
#include <iostream>
#include <locale>

int main() {
    // Установка русской локали
    std::locale::global(std::locale("ru_RU.UTF-8"));
    std::wcout.imbue(std::locale("ru_RU.UTF-8"));

    std::cout << "Решение уравнения Пуассона" << std::endl;
    std::cout << "==========================================================" << std::endl;

    // Параметры сетки
    int N1 = 16, N2 = 16;
    double L1 = 1.0, L2 = 1.0;

    // Создаем решатель
    PoissonSolver2D solver(N1, N2, L1, L2);

    std::cout << "Параметры сетки:" << std::endl;
    std::cout << "  N1 = " << N1 << ", N2 = " << N2 << std::endl;
    std::cout << "  h1 = " << solver.get_hx() << ", h2 = " << solver.get_hy() << std::endl;
    std::cout << "  Общее число точек: " << solver.get_Nx() + 1 << " x " << solver.get_Ny() + 1 << std::endl;

    // Создаем правую часть
    std::vector<std::vector<double>> f(N1 - 1, std::vector<double>(N2 - 1, 0.0));
    for (int i = 0; i < N1 - 1; i++) {
        for (int j = 0; j < N2 - 1; j++) {
            double x = (i + 1) * solver.get_hx();
            double y = (j + 1) * solver.get_hy();
            f[i][j] = PoissonSolver2D::f_test(x, y);
        }
    }

    // Аналитическое решение
    auto analytical = PoissonSolver2D::analytical_solution(N1, N2, L1, L2);

    std::cout << "\nРешаем уравнение Пуассона различными методами:" << std::endl;

    // Решение DST методом
    std::cout << "\n1. Метод дискретного синус-преобразования (DST):" << std::endl;
    auto solution_DST = solver.solve_DST(f);
    double max_error_DST = PoissonSolver2D::compute_max_error(solution_DST, analytical);
    double mean_error_DST = PoissonSolver2D::compute_mean_error(solution_DST, analytical);
    std::cout << "   Максимальная погрешность: " << max_error_DST << std::endl;
    std::cout << "   Средняя погрешность: " << mean_error_DST << std::endl;
    PoissonSolver2D::write_to_file(solution_DST, "solution_DST.txt", L1, L2);
    /*
    // Решение ADI методом
    std::cout << "\n2. Метод переменных направлений (ADI):" << std::endl;
    auto solution_ADI = solver.solve_ADI(f);
    double max_error_ADI = PoissonSolver2D::compute_max_error(solution_ADI, analytical);
    double mean_error_ADI = PoissonSolver2D::compute_mean_error(solution_ADI, analytical);
    std::cout << "   Максимальная погрешность: " << max_error_ADI << std::endl;
    std::cout << "   Средняя погрешность: " << mean_error_ADI << std::endl;
    PoissonSolver2D::write_to_file(solution_ADI, "solution_ADI.txt", L1, L2);
    */
    // Решение попеременно-треугольным методом
    std::cout << "\n3. Попеременно-треугольный метод:" << std::endl;
    auto solution_triangular = solver.solve_triangular_method(f);
    double max_error_triangular = PoissonSolver2D::compute_max_error(solution_triangular, analytical);
    double mean_error_triangular = PoissonSolver2D::compute_mean_error(solution_triangular, analytical);
    std::cout << "   Максимальная погрешность: " << max_error_triangular << std::endl;
    std::cout << "   Средняя погрешность: " << mean_error_triangular << std::endl;
    PoissonSolver2D::write_to_file(solution_triangular, "solution_triangular.txt", L1, L2);

    // Решение шахматным методом
    std::cout << "\n4. Шахматный метод:" << std::endl;
    auto solution_checkerboard = solver.solve_checkerboard_method(f);
    double max_error_checkerboard = PoissonSolver2D::compute_max_error(solution_checkerboard, analytical);
    double mean_error_checkerboard = PoissonSolver2D::compute_mean_error(solution_checkerboard, analytical);
    std::cout << "   Максимальная погрешность: " << max_error_checkerboard << std::endl;
    std::cout << "   Средняя погрешность: " << mean_error_checkerboard << std::endl;
    PoissonSolver2D::write_to_file(solution_checkerboard, "solution_checkerboard.txt", L1, L2);

    // Запись аналитического решения
    PoissonSolver2D::write_to_file(analytical, "analytical_solution.txt", L1, L2);

    std::cout << "\n==========================================================" << std::endl;
    std::cout << "Все решения записаны в файлы:" << std::endl;
    std::cout << "  - solution_DST.txt" << std::endl;
    std::cout << "  - solution_ADI.txt" << std::endl;
    std::cout << "  - solution_triangular.txt" << std::endl;
    std::cout << "  - solution_checkerboard.txt" << std::endl;
    std::cout << "  - analytical_solution.txt" << std::endl;

    return 0;
}