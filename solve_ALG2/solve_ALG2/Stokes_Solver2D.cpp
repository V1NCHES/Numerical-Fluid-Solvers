#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>  
#include "Stokes_Solver2D.h"
#include "Slipper_Area.h"
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <chrono>  


void StokesSolver2D::create_initial_conditions()
{
    std::cout << "Создание матрицы fv1, fv2" << "; Размер f1: " << (nx1 + 1 + nx2 + 1) << " x " << (ny1 + 1) << "Размер f2: " << (nx1 + 1) << " x " << (ny2) << std::endl;

    for (int i = 0; i < nx_d; ++i)
        v2[i][ny2 - 1] = -0.1;

    for (int j = 1; j < ny1; ++j)
        for (int i = 1; i < nx1 + nx2; ++i)
            fu1[i][j] = 1.0, fv1[i][j] = 1.0;

    for (int i = 1; i < nx1; ++i)
        fu1[i][ny1] = 1.0, fv1[i][ny1] = 1.0;

    for (int j = 0; j < ny2 - 1; ++j)
        for (int i = 1; i < nx1; ++i)
            fu2[i][j] = 1.0, fv2[i][j] = 1.0;

    //print_v(f1, f2, nx1 + 1 + nx2, nx1 + 1, ny1, ny2 - 1);
    std::cout << "Матрица f создана" << std::endl;
}



StokesSolver2D::StokesSolver2D(double hx, double hy, double Lx1, double Lx2, double Ly1, double Ly2, double r, double mu, double tau_s, double alpha, const std::string& filename)
{
    this->hx = hx;  this->hy = hy;

    this->Lx1 = Lx1, this->Lx2 = Lx2;  this->Ly1 = Ly1, this->Ly2 = Ly2;
    this->nx1 = Lx1 / hx, this->nx2 = Lx2 / hx;  this->ny1 = Ly1 / hy, this->ny2 = Ly2 / hy;

    this->d_t = 0.001;  this->h2 = hy * hy;
    double sigma = 2 * this->h2 / this->d_t;  this->sm2 = sigma - 2; this->osp2 = 1 / (sigma + 2);

    this->r = r; this->tau_s = tau_s; this->mu = mu; this->alpha = alpha; this->current_iteration = 0; this->filename = filename;
    this->nx_d = nx1 / 2;

    v1 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 1, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2 + 1) x ( ny1 + 1 )
    v2 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)
    u1 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 1, 0.0));       
    u2 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));       

    fv1 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 1, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2 + 1) x ( ny1 + 1 )
    fv2 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)
    fu1 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 1, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2 + 1) x ( ny1 + 1 )
    fu2 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)

    Fv1 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 1, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2 + 1) x ( ny1 + 1 )
    Fv2 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)
    Fu1 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 1, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2 + 1) x ( ny1 + 1 )
    Fu2 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)

    p1 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 1, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2 + 1) x ( ny1 + 1 )
    p2 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)

    g1 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 1, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2 + 1) x ( ny1 + 1 )
    g2 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));       // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)
    r1 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 1, 0.0));
    r2 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));

    this->create_initial_conditions();

    //save_V_file(f1, f2, filename);

    reset_timers();   // Инициализация таймеров

    //grad_v11 = std::vector<std::vector<double>>(nx1 + nx2, std::vector<double>(ny1 + 1, 0.0)); // По x v1
    //grad_v12 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1, 0.0));// По y v1
    //grad_v21 = std::vector<std::vector<double>>(nx1, std::vector<double>(ny2, 0.0));// По x v2
    //grad_v22 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));// По y v2

    //lambda11 = std::vector<std::vector<double>>(nx1 + nx2, std::vector<double>(ny1 + 1, 0.0)); // 
    //lambda12 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1, 0.0));// 
    //lambda21 = std::vector<std::vector<double>>(nx1, std::vector<double>(ny2, 0.0));// 
    //lambda22 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));// 

    //q11 = std::vector<std::vector<double>>(nx1 + nx2, std::vector<double>(ny1 + 1, 0.0)); // 
    //q12 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1, 0.0));// 
    //q21 = std::vector<std::vector<double>>(nx1, std::vector<double>(ny2, 0.0));// 
    //q22 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));// 

    //div_lambda1 = std::vector<std::vector<double>>(nx1 + nx2 - 1, std::vector<double>(ny1 - 1, 0.0));
    //div_lambda2 = std::vector<std::vector<double>>(nx1 - 1, std::vector<double>(ny2, 0.0));
    //div_q1 = std::vector<std::vector<double>>(nx1 + nx2 - 1, std::vector<double>(ny1 - 1, 0.0));
    //div_q2 = std::vector<std::vector<double>>(nx1 - 1, std::vector<double>(ny2, 0.0));

    //t11 = std::vector<std::vector<double>>(nx1 + nx2, std::vector<double>(ny1 + 1, 0.0)); // 
    //t12 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1, 0.0));// 
    //t21 = std::vector<std::vector<double>>(nx1, std::vector<double>(ny2, 0.0));// 
    //t22 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));// 

    //t11_extended = std::vector<std::vector<double>>(nx1 + nx2 + 2, std::vector<double>(ny1 + 1, 0.0)); // 
    //t12_extended = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1 + 2, 0.0));// 
    //t21_extended = std::vector<std::vector<double>>(nx1 + 2, std::vector<double>(ny2, 0.0));// 
    //t22_extended = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2 + 2, 0.0));// 


    //t11_avg = std::vector<std::vector<double>>(nx1 + nx2 + 1, std::vector<double>(ny1, 0.0));
    //t12_avg = std::vector<std::vector<double>>(nx1 + nx2, std::vector<double>(ny1 + 1, 0.0));;
    //t21_avg = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));
    //t22_avg = std::vector<std::vector<double>>(nx1, std::vector<double>(ny2, 0.0));

    //norma11 = std::vector<std::vector<double>>(nx1 + nx2, std::vector<double>(ny1 + 1, 0.0)); // По x v1
    //norma12 = std::vector<std::vector<double>>(nx1 + 1 + nx2, std::vector<double>(ny1, 0.0));// По y v1
    //norma21 = std::vector<std::vector<double>>(nx1, std::vector<double>(ny2, 0.0));// По x v2
    //norma22 = std::vector<std::vector<double>>(nx1 + 1, std::vector<double>(ny2, 0.0));// По y v2


    //norma_lambda11 = std::vector<std::vector<double>>(nx1 + nx2, std::vector<double>(ny1 + 1, 0.0)); // По x v1
    //norma_lambda21 = std::vector<std::vector<double>>(nx1, std::vector<double>(ny2, 0.0));// По x v2

    /*std::ofstream file_diff(this->filename + "diff.txt");
    if (!file_diff.is_open()) { std::cout << "Ошибка открытия файла: " << this->filename + "diff.txt" << std::endl; return; }
    file_diff.close();*/
}


StokesSolver2D::StokesSolver2D()
{
    this->hx = 0.01; this->hy = 0.01; this->Lx1 = 1.0, this->Lx2 = 1.0; this->Ly1 = 1.0, this->Ly2 = 1.0;
}
StokesSolver2D::~StokesSolver2D()
{
    std::cout << " Все почистели\n";
}

double StokesSolver2D::U_ij(double u_ij, double u_ip1j, double u_ijp1, double u_im1j, double u_ijm1, double f_ij)
{
    return (this->sm2 * u_ij + u_ip1j + u_im1j + u_ijp1 + u_ijm1 + f_ij * this->h2) * this->osp2;
}
double StokesSolver2D::Error(double u_ij, double u_ip1j, double u_ijp1, double u_im1j, double u_ijm1, double f_ij)
{
    return fabs(-4.0 * u_ij + u_ip1j + u_im1j + u_ijp1 + u_ijm1 + f_ij * this->h2);
}

double StokesSolver2D::Error_Area(const std::vector<std::vector<double>>& data1, const std::vector<std::vector<double>>& data2)
{
    double Error_max = fabs(data1[1][1] * h2);

    for (int j = 1; j < ny1; ++j)
        for (int i = 1; i < nx1 + nx2; ++i)
            Error_max = std::max(Error_max, fabs(data1[i][j] * h2));

    for (int i = 1; i < nx1; ++i)
        Error_max = std::max(Error_max, fabs(data1[i][ny1] * h2));

    for (int j = 0; j < ny2 - 1; ++j)
        for (int i = 1; i < nx1; ++i)
            Error_max = std::max(Error_max, fabs(data2[i][j] * h2));

    return Error_max;
}



void StokesSolver2D::PoissonSolver2D(double EPS, int max_iter, std::vector<std::vector<double>>& V1, std::vector<std::vector<double>>& V2, std::vector<std::vector<double>>& F1, std::vector<std::vector<double>>& F2)
{
    std::cout << std::endl << "============================================================" << std::endl;
    std::cout << "ШАГ 1: Решение уравнения Пуассона" << std::endl;
    std::cout << "============================================================" << std::endl;

    double Error_max_f = Error_Area(F1, F2);
    std::cout << "Error_max_f = " << Error_max_f << "; EPS * Error_max_f =  " << EPS * Error_max_f << "\n";
    double Error_max_v = 0;
    int k = 0;
    do {
        // n -> n + 1/2
        for (int j = 1; j < ny1; ++j)
            for (int i = 1; i < nx1 + nx2; ++i)
                V1[i][j] = U_ij(V1[i][j], V1[i + 1][j], V1[i][j + 1], V1[i - 1][j], V1[i][j - 1], F1[i][j]);

        for (int i = 1; i < nx1; ++i)
            V1[i][ny1] = U_ij(V1[i][ny1], V1[i + 1][ny1], V2[i][0], V1[i - 1][ny1], V1[i][ny1 - 1], F1[i][ny1]);

        for (int i = 1; i < nx1; ++i)
            V2[i][0] = U_ij(V2[i][0], V2[i + 1][0], V2[i][1], V2[i - 1][0], V1[i][ny1], F2[i][0]);

        for (int j = 1; j < ny2 - 1; ++j)
            for (int i = 1; i < nx1; ++i)
                V2[i][j] = U_ij(V2[i][j], V2[i + 1][j], V2[i][j + 1], V2[i - 1][j], V2[i][j - 1], F2[i][j]);

        // n + 1/2 -> n + 1

        for (int j = ny2 - 2; j > 1; --j)
            for (int i = nx1 - 1; i > 0; --i)
                V2[i][j] = U_ij(V2[i][j], V2[i + 1][j], V2[i][j + 1], V2[i - 1][j], V2[i][j - 1], F2[i][j]);

        for (int i = nx1 - 1; i > 0; --i)
            V2[i][0] = U_ij(V2[i][0], V2[i + 1][0], V2[i][1], V2[i - 1][0], V1[i][ny1], F2[i][0]);

        for (int i = nx1 - 1; i > 0; --i)
            V1[i][ny1] = U_ij(V1[i][ny1], V1[i + 1][ny1], V2[i][0], V1[i - 1][ny1], V1[i][ny1 - 1], F1[i][ny1]);

        for (int j = ny1 - 1; j > 0; --j)
            for (int i = nx1 + nx2 - 1; i > 0; --i)
                V1[i][j] = U_ij(V1[i][j], V1[i + 1][j], V1[i][j + 1], V1[i - 1][j], V1[i][j - 1], F1[i][j]);
  
        Error_max_v = 0.0;

        for (int j = 1; j < ny1; ++j)
            for (int i = 1; i < nx1 + nx2; ++i)
                Error_max_v = std::max(Error_max_v, Error(V1[i][j], V1[i + 1][j], V1[i][j + 1], V1[i - 1][j], V1[i][j - 1], F1[i][j]));

        for (int i = 1; i < nx1; ++i)
            Error_max_v = std::max(Error_max_v, Error(V1[i][ny1], V1[i + 1][ny1], V2[i][0], V1[i - 1][ny1], V1[i][ny1 - 1], F1[i][ny1]));

        for (int i = 1; i < nx1; ++i)
            Error_max_v = std::max(Error_max_v, Error(V2[i][0], V2[i + 1][0], V2[i][1], V2[i - 1][0], V1[i][ny1], F2[i][0]));

        for (int j = 1; j < ny2 - 1; ++j)
            for (int i = 1; i < nx1; ++i)
                Error_max_v = std::max(Error_max_v, Error(V2[i][j], V2[i + 1][j], V2[i][j + 1], V2[i - 1][j], V2[i][j - 1], F2[i][j]));

        if (k % 1000 == 0) std::cout << " k = " << k << "; Error_max_v = " << Error_max_v << "; EPS * Error_max_f = " << EPS * Error_max_f << "\n";
        ++k;
    } while (Error_max_v > EPS * Error_max_f && k < max_iter);
    std::cout << "Final k = " << k << "\n";

    //check_symmetry_detailed(); // ДОБАВИТЬ ПРОВЕРКУ СИММЕТРИИ ЗДЕСЬ
}

void StokesSolver2D::Residual_calculation()
{
    for (int i = 1; i < nx1 + nx2; ++i)
        for (int j = 1; j < ny1; j++)
            g1[i][j] = (u1[i][j] - u1[i-1][j])/hx + (v1[i][j]-v1[i][j-1])/hy;

    for (int i = 1; i < nx1; ++i)
        g1[i][ny1] = (u1[i][ny1] - u1[i - 1][ny1]) / hx + (v2[i][0] - v1[i][ny1]) / hy;

    for (int i = 1; i < nx1; ++i)
        for (int j = 0; j < ny2-1; j++)
            g2[i][j] = (u2[i][j] - u2[i - 1][j]) / hx + (v2[i][j+1] - v2[i][j]) / hy;
}


void  StokesSolver2D::update_Poisson()
{
    //divergence(lambda11, lambda12, div_lambda1, lambda21, lambda22, div_lambda2);
    //divergence(q11, q12, div_q1, q21, q22, div_q2);

    for (int i = 1; i < nx1 + nx2; ++i)
        for (int j = 1; j < ny1; j++)
        {
            Fu1[i][j] = (fu1[i][j]) / r  - (p1[i][j - 1] - p1[i - 1][j - 1]) / hx;
            Fv1[i][j] = (fv1[i][j]) / r  - (p1[i - 1][j] - p1[i - 1][j - 1]) / hy;
        }
            

    for (int i = 1; i < nx1; ++i)
    {
        Fu1[i][ny1] = (fu1[i][ny1]) / r - (p1[i][ny1 - 1] - p1[i - 1][ny1 - 1]) / hx;
        Fv1[i][ny1] = (fv1[i][ny1]) / r - (p2[i - 1][0] - p1[i - 1][ny1]) / hy;
    }


    for (int i = 1; i < nx1; ++i)
        for (int j = 0; j < ny2-1; j++)
        {
            Fu2[i][j] = (fu2[i][j]) / r -(p2[i][j] - p2[i - 1][j]) / hx;
            Fv2[i][j] = (fv2[i][j]) / r - (p2[i - 1][j+1] - p2[i - 1][j]) / hy;
            //std::cout << "i = " << i << "; j = " << j <<std::endl;
        }
           

    //for (int i = 1; i < nx1 + nx2; ++i)
    //    for (int j = 1; j < ny1; j++)
    //        this->F1[i][j] = (f1[i][j] + div_lambda1[i - 1][j - 1] - this->r * this->div_q1[i - 1][j - 1]) / this->r; // f1[i][j] / this-> r;

    //for (int i = 1; i < nx1; ++i)
    //    this->F1[i][ny1] = (f1[i][ny1] + this->div_lambda2[i - 1][0] - this->r * this->div_q2[i - 1][0]) / this->r; //f1[i][ny1] / this->r;

    //for (int i = 1; i < nx1; ++i)
    //    for (int j = 1; j < ny2; j++)
    //        this->F2[i][j] = (f2[i][j] + this->div_lambda2[i - 1][j] - this->r * this->div_q2[i - 1][j]) / this->r;// f2[i][j] / this->r;
}
/*
void StokesSolver2D::divergence(std::vector<std::vector<double>>& field11, std::vector<std::vector<double>>& field12, std::vector<std::vector<double>>& div1, std::vector<std::vector<double>>& field21, std::vector<std::vector<double>>& field22, std::vector<std::vector<double>>& div2)
{
    for (int i = 1; i < nx1 + nx2; ++i)
        for (int j = 1; j < ny1; ++j)
            div1[i - 1][j - 1] = (field11[i][j] - field11[i - 1][j]) / hx + (field12[i][j] - field12[i][j - 1]) / hy;

    //for (int i = 1; i < nx1; ++i)
    //    div1[i - 1][ny1 - 1] = (field11[i][ny1] - field11[i - 1][ny1]) / hx + (field12[i][ny1] - field12[i][ny1 - 1]) / hy;

    for (int i = 1; i < nx1 - 1; ++i)
        div2[i - 1][0] = (field21[i][0] - field21[i - 1][0]) / hx + (field22[i][0] - field12[i][ny1 - 1]) / hy; //?????????????

    for (int i = 1; i < nx1 - 1; ++i)
        for (int j = 1; j < ny2 - 1; ++j)
            div2[i - 1][j] = (field21[i][j] - field21[i - 1][j]) / hx + (field22[i][j] - field22[i][j - 1]) / hy;
}

void StokesSolver2D::gradient_v()
{
    // по x
    for (int i = 0; i < nx1 + nx2; ++i)
        for (int j = 0; j < ny1 + 1; ++j)
            grad_v11[i][j] = (v1[i + 1][j] - v1[i][j]) / hx;

    //for (int i = 0; i < nx1; ++i)
    //    grad_v11[i][ny1] = (v1[i + 1][ny1] - v1[i][ny1]) / hx;

    for (int i = 0; i < nx1; ++i)
        for (int j = 0; j < ny2; ++j)
            grad_v21[i][j] = (v2[i + 1][j] - v2[i][j]) / hx;

    //по y

    for (int i = 0; i < nx1 + nx2 + 1; ++i)
        for (int j = 0; j < ny1; ++j)
            grad_v12[i][j] = (v1[i][j + 1] - v1[i][j]) / hy;

    for (int i = 0; i < nx1 + 1; ++i)
        grad_v22[i][0] = (v2[i][0] - v1[i][ny1]) / hy;

    for (int i = 0; i < nx1 + 1; ++i)
        for (int j = 1; j < ny2; ++j)
            grad_v22[i][j] = (v2[i][j + 1] - v2[i][j]) / hy;

}

void StokesSolver2D::calculation_t(std::vector<std::vector<double>>& lamda, std::vector<std::vector<double>>& grad_v, std::vector<std::vector<double>>& t, std::vector<std::vector<double>>& t_extended, int Nx, int Ny, int k, int l)
{
    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            t[i][j] = lamda[i][j] + this->r * grad_v[i][j], t_extended[i + k][j + l] = t[i][j];
}

void StokesSolver2D::averaging_t(std::vector<std::vector<double>>& t_avg, std::vector<std::vector<double>>& t_extended, int Nx, int Ny)
{
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            t_avg[i][j] = 0.25 * (t_extended[i][j] + t_extended[i + 1][j] + t_extended[i][j + 1] + t_extended[i + 1][j + 1]);
}

void StokesSolver2D::calculation_norma(std::vector<std::vector<double>>& t, std::vector<std::vector<double>>& t_avg, std::vector<std::vector<double>>& norma, int Nx, int Ny)
{
    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
            norma[i][j] = sqrt(t[i][j] * t[i][j] + t_avg[i][j] * t_avg[i][j]);
}

void StokesSolver2D::calculation_q(std::vector<std::vector<double>>& t, std::vector<std::vector<double>>& norma, std::vector<std::vector<double>>& q, int Nx, int Ny)
{
    for (int i = 0; i < Nx; ++i)// {
        for (int j = 0; j < Ny; ++j) //{
            if (norma[i][j] <= this->tau_s) q[i][j] = 0;
            else q[i][j] = t[i][j] * (1 - this->tau_s / norma[i][j]) / (this->mu + this->r);
}

void StokesSolver2D::update_q() {
    std::cout << std::endl << "============================================================" << std::endl;
    std::cout << "ШАГ 2: Обновление q" << std::endl;
    std::cout << "============================================================" << std::endl;

    calculation_t(lambda11, grad_v11, t11, t11_extended, nx1 + nx2, ny1 + 1, 1, 0);
    calculation_t(lambda12, grad_v12, t12, t12_extended, nx1 + 1 + nx2, ny1, 0, 1);
    calculation_t(lambda21, grad_v21, t21, t21_extended, nx1, ny2, 1, 0);
    calculation_t(lambda22, grad_v22, t22, t22_extended, nx1 + 1, ny2, 0, 1);


    averaging_t(t11_avg, t11_extended, nx1 + nx2 + 1, ny1);
    averaging_t(t12_avg, t12_extended, nx1 + nx2, ny1 + 1);
    averaging_t(t21_avg, t21_extended, nx1 + 1, ny2);
    averaging_t(t22_avg, t22_extended, nx1, ny2);


    calculation_norma(t11, t12_avg, norma11, nx1 + nx2, ny1 + 1);
    calculation_norma(t12, t11_avg, norma12, nx1 + 1 + nx2, ny1);
    calculation_norma(t21, t22_avg, norma21, nx1, ny2);
    calculation_norma(t22, t21_avg, norma22, nx1 + 1, ny2);


    calculation_q(t11, norma11, q11, nx1 + nx2, ny1 + 1);
    calculation_q(t12, norma12, q12, nx1 + 1 + nx2, ny1);
    calculation_q(t21, norma21, q21, nx1, ny2);
    calculation_q(t22, norma22, q22, nx1 + 1, ny2);
}

void StokesSolver2D::calculation_lamda(std::vector<std::vector<double>>& lambda, std::vector<std::vector<double>>& grad_v, std::vector<std::vector<double>>& q, int Nx, int Ny)
{
    double lambda1_new = 0.0, lambda2_new = 0.0;
    double diff = 0.0;

    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            lambda1_new = lambda[i][j] + this->alpha * (grad_v[i][j] - q[i][j]);
            diff = fabs(lambda1_new - lambda[i][j]);
            if (diff > this->diff_max) this->diff_max = diff;
            lambda[i][j] = lambda1_new;
            //this->norm_lambda1[i][j] = fabs(this->alpha * (grad_v1[i][j] - this->q1[i][j]));
        }
    }
    //return diff;
}

void StokesSolver2D::update_lambda() {
    std::cout << std::endl << "============================================================" << std::endl;
    std::cout << "ШАГ 3: Обновление lambda" << std::endl;
    std::cout << "============================================================" << std::endl;

    this->diff_max = 0;

    calculation_lamda(lambda11, grad_v11, q11, nx1 + nx2, ny1 + 1);
    calculation_lamda(lambda12, grad_v12, q12, nx1 + 1 + nx2, ny1);
    calculation_lamda(lambda21, grad_v21, q21, nx1, ny2);
    calculation_lamda(lambda22, grad_v22, q22, nx1 + 1, ny2);

    for (int i = 0; i < nx1 + nx2; i++) {
        for (int j = 0; j < ny1 + 1; j++) {
            norma_lambda11[i][j] = fabs(alpha * (grad_v11[i][j] - q11[i][j]));
        }
    }

    for (int i = 0; i < nx1; i++) {
        for (int j = 0; j < ny2; j++) {
            norma_lambda21[i][j] = fabs(alpha * (grad_v21[i][j] - q21[i][j]));
        }
    }

    std::cout << "Максимальное изменение lambda : " << this->diff_max << std::endl;
}

*/

void StokesSolver2D::Two_step_iterative_method()
{
    int k = 0;
    double norma_r1 = 0.0;
    double norma = 0.0;
    while (true)
    {
        for (int i = 1; i < nx1 + nx2; ++i)
            for (int j = 1; j < ny1; j++)
            {
                r1[i][j] = (p1[i][j] - p1[i - 1][j]) / hx + (p1[i][j] - p1[i][j - 1]) / hy - g1[i][j];
                norma = fabs(r1[i][j]);
                if (norma > norma_r1) norma_r1 = norma;
            }               

        for (int i = 1; i < nx1; ++i)
        {
            r1[i][ny1] = (p1[i][ny1] - p1[i - 1][ny1]) / hx + (p2[i][0] - p1[i][ny1]) / hy - g1[i][ny1];
            norma = fabs(r1[i][ny1]);
            if (norma > norma_r1) norma_r1 = norma;
        }
           

        for (int i = 1; i < nx1; ++i)
            for (int j = 0; j < ny2 - 1; j++)
            {
                r2[i][j] = (p2[i][j] - p2[i - 1][j]) / hx + (p2[i][j + 1] - p2[i][j]) / hy - g2[i][j];
                norma = fabs(r2[i][j]);
                if (norma > norma_r1) norma_r1 = norma;
            }

         
        if (norma_r1 < 1e-6) break;
    }
}

int StokesSolver2D::run_full_algorithm(int max_iterations, double convergence_tol) {
    std::cout << "ЗАПУСК ПОЛНОГО АЛГОРИТМА" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "Параметры сетки:" << std::endl;
    std::cout << "  Шаги: hx = " << this->hx << ", hy = " << this->hy << std::endl;
    std::cout << "  Длины по X: Lx1 = " << this->Lx1 << ", Lx2 = " << this->Lx2 << ";  Длины по Y : Ly1 = " << this->Ly1 << ", Ly2 = " << this->Ly2 << std::endl;
    std::cout << "  Количество узлов по X: nx1 = " << this->nx1 << ", nx2 = " << this->nx2 << "  Количество узлов по Y: ny1 = " << this->ny1 << ", ny2 = " << this->ny2 << std::endl;
    std::cout << "Параметры: r = " << this->r << ", mu =" << this->mu << ", tau_s =" << this->tau_s << ", alpha =" << this->alpha << std::endl;
    std::cout << "Максимальное число итераций: " << max_iterations << std::endl;
    std::cout << "Критерий сходимости: " << convergence_tol << std::endl;

    // Переменные для общего времени
    auto total_start_time = std::chrono::high_resolution_clock::now();
    double total_algorithm_time = 0.0;

    int converged = 0;

    for (int iteration = 0; iteration < max_iterations; iteration++) {
        this->current_iteration = iteration + 1;

        std::cout << std::endl << "============================================================" << std::endl;
        std::cout << "ИТЕРАЦИЯ " << this->current_iteration << std::endl;
        std::cout << "============================================================" << std::endl;

        // Сбрасываем таймеры для новой итерации
        reset_timers();
        auto iteration_start = std::chrono::high_resolution_clock::now();

        // Шаг 1: Решение уравнения Пуассона с замером времени
        auto poisson_start = std::chrono::high_resolution_clock::now();
        update_Poisson();

        PoissonSolver2D(1e-6, 10000, u1, u2, Fu1, Fu2);
        PoissonSolver2D(1e-6, 10000, v1, v2, Fv1, Fv2);

        Residual_calculation();

        std::cout << "g = \n";
        print_v(g1, g2, nx1 + 1 + nx2, nx1 + 1, ny1, ny2 - 1);
        std::cout << "\n";

        std::cout << "Fu = \n";
        print_v(Fu1, Fu2, nx1 + 1 + nx2, nx1 + 1, ny1, ny2 - 1);
        std::cout << "\n";

        std::cout << "Fv = \n";
        print_v(Fv1, Fv2, nx1 + 1 + nx2, nx1 + 1, ny1, ny2 - 1);
        /*
        std::cout << "u = \n";
        print_v(u1, u2, nx1 + 1 + nx2, nx1 + 1, ny1, ny2 - 1);
        std::cout << "\n";

        std::cout << "v = \n";
        print_v(v1, v2, nx1 + 1 + nx2, nx1 + 1, ny1, ny2 - 1);
        std::cout << "\n";

        std::cout << "fu = \n";
        print_v(fu1, fu2, nx1 + 1 + nx2, nx1 + 1, ny1, ny2 - 1);
        std::cout << "\n";

        std::cout << "fv = \n";
        print_v(fv1, fv2, nx1 + 1 + nx2, nx1 + 1, ny1, ny2 - 1);
        std::cout << "\n";*/

        
       //PoissonSolver2D(1e-6, 10000);
        auto poisson_end = std::chrono::high_resolution_clock::now();
        poisson_time = std::chrono::duration<double>(poisson_end - poisson_start).count();


        // Шаг 2: Обновление q с замером времени
        auto update_q_start = std::chrono::high_resolution_clock::now();
       /* gradient_v();
        update_q();*/
        auto update_q_end = std::chrono::high_resolution_clock::now();
        update_q_time = std::chrono::duration<double>(update_q_end - update_q_start).count();

        // Шаг 3: Обновление lambda с замером времени
        auto update_lambda_start = std::chrono::high_resolution_clock::now();
        /*update_lambda();*/
        auto update_lambda_end = std::chrono::high_resolution_clock::now();
        update_lambda_time = std::chrono::duration<double>(update_lambda_end - update_lambda_start).count();

        // Общее время итерации
        auto iteration_end = std::chrono::high_resolution_clock::now();
        total_time = std::chrono::duration<double>(iteration_end - iteration_start).count();

        print_iteration_time(); // Выводим время выполнения итерации

        Save();

        if (this->diff_max < convergence_tol) {
            std::cout << std::endl << "СХОДИМОСТЬ ДОСТИГНУТА на итерации " << this->current_iteration << std::endl;
            converged = 1;
            break;
        }

        if (this->current_iteration % 10 == 0) {
            std::cout << std::endl << "--- Прогресс: итерация " << this->current_iteration
                << ", изменение lambda = " << std::scientific << this->diff_max << std::endl;
        }

    }

    // Вычисляем общее время работы алгоритма
    auto total_end_time = std::chrono::high_resolution_clock::now();
    total_algorithm_time = std::chrono::duration<double>(total_end_time - total_start_time).count();

    // Сохраняем финальные данные
    //Save();
    save_V_file(v1, v2, nx1 + 1 + nx2, ny1 + 1, nx1 + 1, ny2, "data1/v_matrix_final.txt");

    // Выводим итоговую статистику по времени
    std::cout << std::endl << "============================================================" << std::endl;
    std::cout << "ИТОГИ ВРЕМЕНИ ВЫПОЛНЕНИЯ:" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "Общее время алгоритма: " << std::fixed << std::setprecision(2) << total_algorithm_time << " сек" << std::endl;
    std::cout << "Количество итераций: " << this->current_iteration << std::endl;
    if (this->current_iteration > 0) {
        std::cout << "Среднее время итерации: " << std::setprecision(4)
            << (total_algorithm_time / this->current_iteration) << " сек" << std::endl;
    }

    std::cout << std::endl << "============================================================" << std::endl;
    std::cout << "АЛГОРИТМ ЗАВЕРШЕН" << std::endl;
    if (converged) {
        std::cout << "АЛГОРИТМ СОШЕЛСЯ за " << this->current_iteration << " итераций" << std::endl;
    }
    else {
        std::cout << "ДОСТИГНУТО МАКСИМАЛЬНОЕ ЧИСЛО ИТЕРАЦИЙ (" << max_iterations << ")" << std::endl;
    }

    return converged;
}

/*


*/


void StokesSolver2D::Save()
{
    if (this->current_iteration <= 200)
    {
        std::string file_u = this->filename + "u_iteration_" + std::to_string(this->current_iteration) + ".txt";
        save_V_file(u1, u2, nx1 + 1 + nx2, ny1 + 1, nx1 + 1, ny2, file_u.c_str());

        std::string file_v = this->filename + "v_iteration_" + std::to_string(this->current_iteration) + ".txt";
        save_V_file(v1, v2, nx1 + 1 + nx2, ny1 + 1, nx1 + 1, ny2, file_v.c_str());


        //std::string file_q1 = this->filename + "q1_iteration_" + std::to_string(this->current_iteration) + ".txt";
        //save_V_file(q11, q21, nx1 + nx2, ny1 + 1, nx1, ny2, file_q1.c_str());

        //std::string file_lambda1 = this->filename + "lambda1_iteration_" + std::to_string(this->current_iteration) + ".txt";
        //save_V_file(lambda11, lambda21, nx1 + nx2, ny1 + 1, nx1, ny2, file_lambda1.c_str());

        //std::string file_norma1 = this->filename + "norma1_iteration_" + std::to_string(this->current_iteration) + ".txt";
        //save_V_file(norma11, norma21, nx1 + nx2, ny1 + 1, nx1, ny2, file_norma1.c_str());

        //std::string file_norm_lambda1 = this->filename + "norma_lambda1_iteration_" + std::to_string(this->current_iteration) + ".txt";
        //save_V_file(norma_lambda11, norma_lambda21, nx1 + nx2, ny1 + 1, nx1, ny2, file_norm_lambda1.c_str());//*/

        /*std::ofstream file_diff(this->filename + "diff.txt", std::ios::app);
        if (!file_diff.is_open()) { std::cout << "Ошибка открытия файла: " << this->filename + "diff.txt" << std::endl; return; }
        file_diff << this->diff_max << std::endl;
        file_diff.close();*/
    }
}

void StokesSolver2D::save_V_file(const std::vector<std::vector<double>>& data1, const std::vector<std::vector<double>>& data2, int Nx1, int Ny1, int Nx2, int Ny2, const std::string& filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл " << filename << std::endl;
        return;
    }
    double x = 0.0, y = 0.0;

    file << std::scientific << std::setprecision(6);

    for (int i = 0; i < Nx1; ++i)
    {
        for (int j = 0; j < Ny1; ++j)
        {
            x = (i)*this->hx; y = (j)*this->hy;
            file << x << " " << y << " " << data1[i][j] << std::endl;
        }
        file << std::endl;
    }
    for (int i = 0; i < Nx2; ++i)
    {
        for (int j = 0; j < Ny2; ++j)
        {
            x = (i)*this->hx; y = (j + Ny1) * this->hy;
            file << x << " " << y << " " << data2[i][j] << std::endl;
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Данные записаны в файл: " << filename << std::endl;
}


//void print_v(std::vector<std::vector<double>>& print1, std::vector<std::vector<double>>& print2, int nx1, int nx2, int ny1, int ny2)
//{
//    for (int j = ny2; j >= 0; --j)
//    {
//        for (int i = 0; i < nx2; ++i)
//        {
//            std::cout << print2[i][j] << " ";
//        }
//        std::cout << "\n";
//    }
//
//    for (int j = ny1; j >= 0; --j)
//    {
//        for (int i = 0; i < nx1; ++i)
//        {
//            std::cout << print1[i][j] << " ";
//        }
//        std::cout << "\n";
//    }
//}


void StokesSolver2D::check_symmetry_detailed() {
    std::cout << "ПРОВЕРКА СИММЕТРИИ:" << std::endl;

    bool symmetric = true;
    double max_diff = 0.0;
    /*
    for (int i = 0; i < this->nx; i++) {
        for (int j = i + 1; j < this->ny; j++) {
            double diff = std::abs(this->v[i][j] - this->v[j][i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
            if (diff > 1e-12) {
                symmetric = false;
            }
        }
    }
    */
    std::cout << "Нужно добавить Матрица v симметрична: " << (symmetric ? "ДА" : "НЕТ") << std::endl;
    std::cout << "Максимальное отклонение от симметрии: " << std::scientific << max_diff << std::endl;
}

// Метод вывода времени итерации
void StokesSolver2D::print_iteration_time() {
    std::cout << "ВРЕМЯ ВЫПОЛНЕНИЯ ИТЕРАЦИИ:" << std::endl;
    std::cout << "  Общее время: " << std::fixed << std::setprecision(4) << total_time << " сек" << std::endl;
    std::cout << "  Решение Пуассона: " << poisson_time << " сек ("
        << std::setprecision(1) << (poisson_time / total_time * 100) << "%)" << std::endl;
    std::cout << "  Обновление q: " << std::setprecision(4) << update_q_time << " сек ("
        << std::setprecision(1) << (update_q_time / total_time * 100) << "%)" << std::endl;
    std::cout << "  Обновление lambda: " << std::setprecision(4) << update_lambda_time << " сек ("
        << std::setprecision(1) << (update_lambda_time / total_time * 100) << "%)" << std::endl;
}