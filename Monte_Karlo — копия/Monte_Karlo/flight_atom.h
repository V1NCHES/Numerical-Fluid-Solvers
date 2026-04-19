#ifndef flight_atom_h
#define flight_atom_h


#include"include.h"
#include <omp.h>
#include"velocity_at_border.h"
#include "particle_collision.h"
#include "free_path_length.h"
#include"distribution_function.h"
using namespace std;
void flight_atom(ÑMatrix& M, int& ptr, int N, vector<vector<double>> &graficV);
void distagramma_normal_extractionem();
void distagramma_velociti_at_border();
void distagramma_free_path_length_hard_balls();
void distagramma_free_path_length_Reloads();
void distagramma_velociti_proton();
void distagramma_new_velociti();
double G(double wx, double wy, double wz, double vx, double vy, double vz);
double f_P(double wx, double wy, double wz, double vx, double vy, double vz, double U);

double method_of_rectangles_R3(double(*f)(double, double, double, double, double, double, double), vector <double> min_lim, vector <double> max_lim, double delta, double v_x, double v_y, double v_z);
double integrate_R3(double(*f)(double, double, double, double, double, double, double), vector <double> min_lim, vector <double> max_lim, int n, double v_x, double v_y, double v_z);
#endif