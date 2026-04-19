#include "LShapeSolver2D.h"
#include <iostream>
#include <windows.h>
#include <string>

int main() {
    // Directories to create
    std::string output = "output_data_boost";
    std::string SAVE = "SAVE ITERATION";
    std::string SAVE_boost = "SAVE ITERATION Boost";
    
    // Create the output directory if it doesn't exist
    CreateDirectoryA(output.c_str(), NULL);
    CreateDirectoryA(SAVE.c_str(), NULL);
    CreateDirectoryA(SAVE_boost.c_str(), NULL);

    LShapeConfig cfg;
    cfg.h = 0.01;
    cfg.nx1 = 100;
    cfg.nx2 = 200;
    cfg.ny1 = 100;
    cfg.ny2 = 200;
    cfg.v_inflow = -0.1;
    cfg.lx_hole = 50; // Half of the top width is the hole

    //std::cout << "Initializing L-Shape Solver..." << std::endl;
    LShapeSolver2D solver(cfg);
    
    // Bingham properties
    solver.mu = 1.0;
    solver.r_param = 1.0;
    solver.rho = 1.0;
    solver.tau_s = 0.1 * std::sqrt(2.0); // Example initial value from Caverna
    
    int max_iter = 50;
    double tol = 1e-4;

    std::cout << "Running full algorithm..." << std::endl;
    solver.run_full_algorithm(max_iter, tol);

    std::cout << "Done." << std::endl;
    return 0;
}
