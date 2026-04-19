#include <iostream>
#include "algorithm_checker.h"
#include "Slipper_Area.h"
#include "Stokes_Solver2D.h"
#include <iostream>
#include <locale>
#include <windows.h>

void main_Algoritm_rectangle()
{
    setlocale(LC_ALL, "Russian");
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
    std::cout << "Запуск алгоритма проверки..." << std::endl;
    const std::string filename = "data_rectangle/data000/";

    AlgorithmChecker checker(0.005, 0.005, 2.0, 1.0, 1.0, 1.0, 0.2, 1.0, filename);
    int result = checker.run_full_algorithm(100, 1e-6);

    if (result) {
        std::cout << "Программа завершена успешно (сходимость достигнута)" << std::endl;
    }
    else {
        std::cout << "Программа завершена (достигнуто максимальное число итераций)" << std::endl;
    }
}

void main_Algoritm_square()
{
    setlocale(LC_ALL, "Russian");
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
    std::cout << "Запуск алгоритма проверки..." << std::endl;
    const std::string filename = "data_square/data2/";

    AlgorithmChecker checker(0.005, 0.005, 1.0, 1.0, 1.0, 1.0, 0.2, 1.0, filename);
    int result = checker.run_full_algorithm(500, 1e-6);

    if (result) {
        std::cout << "Программа завершена успешно (сходимость достигнута)" << std::endl;
    }
    else {
        std::cout << "Программа завершена (достигнуто максимальное число итераций)" << std::endl;
    }
    checker.Save(0);
}

void main_ALG2_slipper()
{
    setlocale(LC_ALL, "Russian"); // Установка русской локали для консоли 
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
    std::cout << "Запуск алгоритма проверки..." << std::endl;
    const std::string filename = "data_L_aera/data000/";
    double hx = 0.005, hy = 0.005;
    double Lx1 = 1.0, Lx2 = 1.0, Ly1 = 1.0, Ly2 = 1.0;

    ALG2_Slipper checker(hx, hy, Lx1, Lx2, Ly1, Ly2, 1.0, 1.0, 0.1, 1.0, filename);

    int result = checker.run_full_algorithm(500, 1e-6);

    if (result) {
        std::cout << "Программа завершена успешно (сходимость достигнута)" << std::endl;
    }
    else {
        std::cout << "Программа завершена (достигнуто максимальное число итераций)" << std::endl;
    }//*/
}

void main_ALG2_Stokes()
{
    setlocale(LC_ALL, "Russian"); // Установка русской локали для консоли 
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
    std::cout << "Запуск алгоритма проверки..." << std::endl;
    const std::string filename = "data3/";
    //double hx = 0.005, hy = 0.005;
    double hx = 0.25, hy = 0.25;
    double Lx1 = 1.0, Lx2 = 1.0, Ly1 = 1.0, Ly2 = 1.0;

    StokesSolver2D checker(hx, hy, Lx1, Lx2, Ly1, Ly2, 1.0, 1.0, 0.2, 1.0, filename);

    int result = checker.run_full_algorithm(1, 1e-6);

    if (result) {
        std::cout << "Программа завершена успешно (сходимость достигнута)" << std::endl;
    }
    else {
        std::cout << "Программа завершена (достигнуто максимальное число итераций)" << std::endl;
    }//*/
}


int main() {
    
    //main_Algoritm_square();

   // main_Algoritm_rectangle();

   main_ALG2_slipper();
    //main_ALG2_Stokes();

    return 0;
}