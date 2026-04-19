#ifndef include_h
#define include_h

#define Stream = omp_set_num_threads(omp_get_num_procs());

#include <windows.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <random>
#include <ctime> 
#include <math.h> 
#include <limits>

using namespace std;


const double inf = numeric_limits<double>::infinity();;
const double PI = acos(-1.0); // число pi
const int N_ = 1e10;
// Попробуй с_Н=1, с_р=1.75098, u_p=-0.963449, u_H=-2.54351
//const double k = 1.380649 * 1e-23 ; //км/с
//const double Tp = 6000; //K
const double k = 1.38064852 * 1e-16;//  протона в граммах
const double mp = 1.6735575e-24; // г
const double T_H = 6530; //K
const double k_mh = 13806.49 / 1.6735575;
const double a1 = 2.2835 * 1e-7;
const double a2 = 1.062 * 1e-8;
const double C_H = sqrt(2* k*T_H/mp);
const double E = 1.6 * 1e-12;
const double Ve = sqrt(4*E/mp);
const double sigma_Ve = fabs(a1 - a2 * log(Ve));
//const double np = 0.04; //cм^(-3)
//const double mp = 1.67262192 * 1e-27;//кг

const double L = 5;//
const double Cp = 1.75; //1.75
const double CH = 1.0;
const double U_H = -1.58;  // -1.58; //км/с  -2.54
const double U_p = 0.0; //км/с  -0.96
//const double l_dimensionless = 1 / (np);
//вот эти нашел из файла
double sigma_(double v);

double extractionem_0_1();
double extractionem_0_1_();
double normal_extractionem_0_1(int &a1_, int &a2_, int &a3_);
double extractionem_1_1();
double getCPUTime();



bool is_equal(double x, double y);
double sqrt_check(double x, int index);
bool check_inf(double V);
double log_check(double x);
#endif