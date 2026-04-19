#include <iostream>
#include <windows.h>
#include <locale>
#include "CavernaSolver2D.h"



int main() {
    // Windows console encoding for Cyrillic support
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    setlocale(LC_ALL, "Russian");

    std::cout << "\nStarting 2D Caverna Solver (Consolidated)...\n" << std::endl;

    double h = 0.005;      // Grid step
    double Lx = 1.0;     // Cavity width
    double Ly = 1.0;     // Cavity height
    double mu = 1.0;     // Viscosity
    double tau_s = sqrt(2);  // Yield stress
    double r = 1.0;      // Algorithm parameter r
    double rho = 1.0;    // Algorithm parameter rho
    std::string output = "output_data_boost";
    std::string SAVE = "SAVE ITERATION"; 
    std::string SAVE_boost = "SAVE ITERATION Boost";
    // Create the output directory if it doesn't exist
    CreateDirectoryA(output.c_str(), NULL);
    CreateDirectoryA(SAVE.c_str(), NULL);
    CreateDirectoryA(SAVE_boost.c_str(), NULL);

    CavernaSolver2D solver(h, Lx, Ly, mu, tau_s, r, rho, output);

    
    //int status = solver.run_full_algorithm(400, 1e-6);
    // 
    int status = solver.run_full_algorithm_boost(400, 1e-6);
    //solver.boundary_conditions();
    /*Timer timer_cg_custom;
    solver.solve_pressure_cg_custom();
    double t_sg = timer_cg_custom.elapsed();*/
    //solver.save_all(1);


    return 0;
}
