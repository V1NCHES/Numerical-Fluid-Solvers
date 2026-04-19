#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>  // Добавить для sprintf
#include "algorithm_checker.h"
#include "DST.h"
#include "PoissonSolver2D.h"
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <chrono>  // Добавить в includes
#include <iostream>
#include <locale>
#include <windows.h>

// Метод сброса таймеров
void AlgorithmChecker::reset_timers() {
    total_time = 0.0;
    poisson_time = 0.0;
    update_q_time = 0.0;
    update_lambda_time = 0.0;
}

// Метод вывода времени итерации
void AlgorithmChecker::print_iteration_time() {
    std::cout << "ВРЕМЯ ВЫПОЛНЕНИЯ ИТЕРАЦИИ:" << std::endl;
    std::cout << "  Общее время: " << std::fixed << std::setprecision(4) << total_time << " сек" << std::endl;
    std::cout << "  Решение Пуассона: " << poisson_time << " сек ("
        << std::setprecision(1) << (poisson_time / total_time * 100) << "%)" << std::endl;
    std::cout << "  Обновление q: " << std::setprecision(4) << update_q_time << " сек ("
        << std::setprecision(1) << (update_q_time / total_time * 100) << "%)" << std::endl;
    std::cout << "  Обновление lambda: " << std::setprecision(4) << update_lambda_time << " сек ("
        << std::setprecision(1) << (update_lambda_time / total_time * 100) << "%)" << std::endl;
}
//Находим максимум матрицы
double max_matrix(double** matrix, int nx, int ny)
{
    double max = matrix[0][0];

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0;j < ny; ++j)
        {
            if (matrix[i][j] > max)
            {
                max = matrix[i][j];
            }
        }
    }
    return max;
}


void AlgorithmChecker::check_symmetry_detailed() {
    std::cout << "ПРОВЕРКА СИММЕТРИИ:" << std::endl;

    bool symmetric = true;
    double max_diff = 0.0;

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

    std::cout << "Матрица v симметрична: " << (symmetric ? "ДА" : "НЕТ") << std::endl;
    std::cout << "Максимальное отклонение от симметрии: " << std::scientific << max_diff << std::endl;
    
}

void AlgorithmChecker::create_f_matrix() {
    // Создает матрицу f: внутренние точки = 1, граничные 0
    // Размер: (nx-2) x (ny-2) - только внутренние точки

    std::cout << "Создание матрицы f..." << std::endl;
    std::cout << "Размер f: " << (nx - 2) << " x " << (ny - 2) << std::endl;

    // Матрица f уже выделена в конструкторе, просто заполняем значениями
    for (int i = 0; i < nx - 2; i++) {
        for (int j = 0; j < ny - 2; j++) {
            f[i][j] = 1.0;  // Все внутренние точки = 1.0
        }
    }

    std::cout << "Матрица f создана (все значения = 1.0)" << std::endl;
}

AlgorithmChecker::AlgorithmChecker(double hx, double hy, double Lx, double Ly, double r, double mu, double tau_s, double alpha, const std::string& filename) {
    this->hx =hx;
    this->hy =hy;
    this->Lx = Lx;
    this->Ly = Ly;
    this->nx = this->Lx / hx  + 1;
    this->ny = this->Ly / hy + 1;

    // Инициализация матриц нулями
    this->lambda1 = allocate_matrix(nx - 1, ny);
    this->lambda2 = allocate_matrix(nx, ny - 1);
    this->q1 = allocate_matrix(nx - 1, ny);
    this->q2 = allocate_matrix(nx, ny - 1);
    this->v = allocate_matrix(nx, ny);
    this->v_f = allocate_matrix(nx - 2, ny - 2);
    
    this-> f_i = allocate_matrix(nx - 2, ny - 2);
    this-> div_lambda = allocate_matrix(nx - 2, ny - 2);
    this-> div_q = allocate_matrix(nx - 2, ny - 2);

    this->F_i = allocate_matrix(nx , ny);

    //
    this->grad_v1 = allocate_matrix(this->nx - 1, this->ny);
    this->grad_v2 = allocate_matrix(this->nx, this->ny - 1);
    this->t1 = allocate_matrix(this->nx - 1, this->ny);
    this->t2 = allocate_matrix(this->nx, this->ny - 1);
    this->t1_extended = allocate_matrix(this->nx + 1, this->ny);
    this->t2_extended = allocate_matrix(this->nx, this->ny + 1);
    this->t1_avg = allocate_matrix(this->nx, this->ny - 1);
    this->t2_avg = allocate_matrix(this->nx - 1, this->ny);
    this->norm1 = allocate_matrix(this->nx - 1, this->ny);
    this-> norm2 = allocate_matrix(this->nx, this->ny - 1);


   // this->norm_lambda1 = allocate_matrix(nx - 1, ny);


    this->r = r;
    this->tau_s = tau_s;
    this->mu = mu;
    this->alpha = alpha;
    this->current_iteration = 0;

    this->filename = filename;
    // Тестовые данные f
    this->f = allocate_matrix(nx - 2, ny - 2);
    this->create_f_matrix();
   

    // Инициализируем решатель
    //this->poisson_solver = new CorrectPoissonSolver2D(nx - 1, ny - 1, this->Lx, this->Ly);
    this->poisson_solver = new PoissonSolver2D(nx - 1, ny - 1, this->Lx, this->Ly);
    
    // Инициализация таймеров
    reset_timers();

    // Координаты сетки
    this->x = new double[nx];
    this->y = new double[ny];
    for (int i = 0; i < nx; i++) {
        this->x[i] = i * this->hx;
    }
    for (int j = 0; j < ny; j++) {
        this->y[j] = j * this->hy;
    }

    // Инициализация истории
    this->convergence_history = nullptr;
    this->v_history = nullptr;
    this->history_size = 0;
}

AlgorithmChecker::~AlgorithmChecker() {
    free_matrix(this->lambda1, this->nx - 1);
    free_matrix(this->lambda2, this->nx);
    free_matrix(this->q1, this->nx - 1);
    free_matrix(this->q2, this->nx);
    free_matrix(this->v, this->nx);
    free_matrix(this->f, this->nx - 2);
    free_matrix(this->v_f, this->nx - 2);
   
    free_matrix(this->f_i, this->nx - 2);
    free_matrix(this->div_lambda, this->nx - 2);
    free_matrix(this->div_q, this->nx - 2);
    free_matrix(this->F_i, this->nx);


    free_matrix(grad_v1, this->nx - 1);
    free_matrix(grad_v2, this->nx);
    free_matrix(t1, this->nx - 1);
    free_matrix(t2, this->nx);
    free_matrix(t1_extended, this->nx + 1);
    free_matrix(t2_extended, this->nx);
    free_matrix(t1_avg, this->nx);
    free_matrix(t2_avg, this->nx - 1);
    free_matrix(norm1, this->nx - 1);
    free_matrix(norm2, this->nx);

    //free_matrix(this->norm_lambda1,this->nx - 1);

    delete[] this->x;
    delete[] this->y;
    //delete this->poisson_solver;
    delete this->poisson_solver;

    if (this->convergence_history) delete[] this->convergence_history;

    if (this->v_history) {
        for (int i = 0; i < this->history_size; i++) {
            free_matrix(this->v_history[i], this->nx);
        }
        delete[] this->v_history;
    }
}

double** AlgorithmChecker::allocate_matrix(int rows, int cols) {
    double** matrix = new double* [rows];
    for (int i = 0; i < rows; i++) {
        matrix[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = 0.0;
        }
    }
    return matrix;
}

void AlgorithmChecker::free_matrix(double** matrix, int rows) {
    for (int i = 0; i < rows; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

void AlgorithmChecker::save_matrix_to_file(double** matrix, int rows, int cols, const char* filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cout << "Ошибка открытия файла: " << filename << std::endl;
        return;
    }
    double x = 0.0, y = 0.0;
    file.precision(6);
    file << std::fixed;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            x = (i)*this->hx; y = (j)*this->hy;
            file << x << " " << y << " " << matrix[i][j] << std::endl;
        }
        file << std::endl;
    }

    file.close();
    std::cout << "Матрица сохранена: " << filename << std::endl;
}

void AlgorithmChecker::save_matrix_to_file_norma(double** norna1, double** norma2, int rows, int cols, const char* filename) 
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cout << "Ошибка открытия файла: " << filename << std::endl;
        return;
    }
    double x = 0.0, y = 0.0;
    file.precision(6);
    file << std::fixed;

    for (int j = 1; j <= ny - 1; j++)  
    {

        for (int i = 1; i <= nx - 1; i++)
        {
            x = (i)*this->hx - hx/2.0; y = (j-1)*this->hy;
            file << x << " " << y << " " << norna1[i-1][j-1] << std::endl;
        }
        file << std::endl;

        for (int i = 0; i <= nx - 1; i++)
        {
            x = (i)*this->hx; y = (j) * this->hy - hy/2;
            file << x << " " << y << " " << norma2[i][j-1] << std::endl;
        }
    }

    for (int i = 1; i <= nx - 1; i++)
    {
        x = (i)*this->hx - hx / 2.0; y = Ly;
        file << x << " " << y << " " << norna1[i-1][ny] << std::endl;
    }

    file.close();
    std::cout << "Матрица сохранена: " << filename << std::endl;
}

void AlgorithmChecker::divergence(double** field1, double** field2, double** div) {
    for (int i = 1; i < this->nx - 1; i++) {
        for (int j = 1; j < this->ny - 1; j++) {
            double div_x = (field1[i][j] - field1[i - 1][j]) / this->hx;
            double div_y = (field2[i][j] - field2[i][j - 1]) / this->hy;
            div[i - 1][j - 1] = div_x + div_y;
        }
    }
}

void AlgorithmChecker::gradient(double** v, double** grad_v1, double** grad_v2) {
    for (int i = 0; i < this->nx - 1; i++) {
        for (int j = 0; j < this->ny; j++) {
            grad_v1[i][j] = (v[i + 1][j] - v[i][j]) / this->hx;
        }
    }

    for (int i = 0; i < this->nx; i++) {
        for (int j = 0; j < this->ny - 1; j++) {
            grad_v2[i][j] = (v[i][j + 1] - v[i][j]) / this->hy;
        }
    }
}

int AlgorithmChecker::solve_poisson_step(double** f) {
    // В функции solve_poisson_step:
    
    std::cout << std::endl << "============================================================" << std::endl;
    std::cout << "ШАГ 1: Решение уравнения Пуассона" << std::endl;
    std::cout << "============================================================" << std::endl;

    divergence(this->lambda1, this->lambda2, this->div_lambda);
    divergence(this->q1, this->q2, this->div_q);

    for (int i = 0; i < this->nx - 2; i++) {
        for (int j = 0; j < this->ny - 2; j++) 
        {
            this->f_i[i][j] = (f[i][j] + this->div_lambda[i][j] - this->r * this->div_q[i][j]) / this->r;
            this->F_i[i + 1][j + 1] = this->f_i[i][j];
        }
    }
    // Решаем уравнение Пуассона 
   // this->poisson_solver->solve(f_i, v_f);
   // this->poisson_solver->solve_DST(this->f_i, this->v_f);
    this->poisson_solver->solve_triangular_method(this->F_i, this->v);
   
    // Заполняем матрицу v
    
    // Заполняем внутренние точки
    /*
    for (int i = 1; i < this->nx - 1; i++) {
        for (int j = 1; j < this->ny - 1; j++) {
            this->v[i][j] = this->v_f[i - 1][j - 1];
        }
    }//*/

    // Устанавливаем граничные условия
    //for (int j = 0; j < this->ny; j++) {
    //    this->v[0][j] = 0.0;      // нижняя граница
    //    this->v[this->nx - 1][j] = 0.0; // верхняя граница
    //}
    //for (int i = 0; i < this->nx; i++) {
    //    this->v[i][0] = 0.0;      // левая граница
    //    this->v[i][this->ny - 1] = 0.0; // правая граница
    //}

    // ДОБАВИТЬ ПРОВЕРКУ СИММЕТРИИ ЗДЕСЬ
    //check_symmetry_detailed();

    // Сохраняем матрицу v каждые 10 итераций
    /*
    if ( this->current_iteration <= 100) { //this->current_iteration % 10 == 0 ||

        std::string file_v = this->filename + "v_iteration_" + std::to_string(this->current_iteration) + ".txt";
        save_matrix_to_file(this->v, this->nx, this->ny, file_v.c_str());

        std::string file_q1 = this->filename + "q1_iteration_" + std::to_string(this->current_iteration) + ".txt";
        save_matrix_to_file(this->q1, this->nx - 1, this->ny, file_q1.c_str());

        std::string file_lambda1 = this->filename + "lambda1_iteration_" + std::to_string(this->current_iteration) + ".txt";
        save_matrix_to_file(this->lambda1, this->nx - 1, this->ny, file_lambda1.c_str());
    }*/
    return 1;
}

void AlgorithmChecker::calculation_norma(double** Q1, double** Q2)
{
    for (int i = 0; i < this->nx - 1; i++) {
        for (int j = 0; j < this->ny; j++) {
            t1_extended[i + 1][j] = Q1[i][j];
        }
    }

    for (int i = 0; i < this->nx; i++) {
        for (int j = 0; j < this->ny - 1; j++) {
            t2_extended[i][j + 1] = Q2[i][j];
        }
    }

    // Упрощенный расчет средних - исправление

    for (int i = 0; i < this->nx; i++) {
        for (int j = 0; j < this->ny - 1; j++) {


            t1_avg[i][j] = 0.25 * (t1_extended[i][j] + t1_extended[i + 1][j] + t1_extended[i][j + 1] + t1_extended[i + 1][j + 1]);
        }
    }

    for (int i = 0; i < this->nx - 1; i++) {
        for (int j = 0; j < this->ny; j++) {

            t2_avg[i][j] = 0.25 * (t2_extended[i][j] + t2_extended[i + 1][j] + t2_extended[i][j + 1] + t2_extended[i + 1][j + 1]);


        }
    }

    // Нормы

    for (int i = 0; i < this->nx - 1; i++) {
        for (int j = 0; j < this->ny; j++) {
            norm1[i][j] = sqrt(Q1[i][j] * Q1[i][j] + t2_avg[i][j] * t2_avg[i][j]);
        }
    }

    for (int i = 0; i < this->nx; i++) {
        for (int j = 0; j < this->ny - 1; j++) {
            norm2[i][j] = sqrt(t1_avg[i][j] * t1_avg[i][j] + Q2[i][j] * Q2[i][j]);
        }
    }
}

void AlgorithmChecker::update_q() {
    std::cout << std::endl << "============================================================" << std::endl;
    std::cout << "ШАГ 2: Обновление q" << std::endl;
    std::cout << "============================================================" << std::endl;

    gradient(this->v, grad_v1, grad_v2);

    for (int i = 0; i < this->nx - 1; i++) {
        for (int j = 0; j < this->ny; j++) {
            t1[i][j] = this->lambda1[i][j] + this->r * grad_v1[i][j];
        }
    }

    for (int i = 0; i < this->nx; i++) {
        for (int j = 0; j < this->ny - 1; j++) {
            t2[i][j] = this->lambda2[i][j] + this->r * grad_v2[i][j];
        }
    }
   
   
    calculation_norma(t1, t2);
   

    // Обновление q
    for (int i = 0; i < this->nx - 1; i++) {
        for (int j = 0; j < this->ny; j++) {
            if (norm1[i][j] <= this->tau_s) {
                this->q1[i][j] = 0;
            }
            else {
                this->q1[i][j] = t1[i][j] * (1 - this->tau_s / norm1[i][j]) / (this->mu + this->r);
            }
        }
    }

    for (int i = 0; i < this->nx; i++) {
        for (int j = 0; j < this->ny - 1; j++) {
            if (norm2[i][j] <= this->tau_s) {
                this->q2[i][j] = 0;
            }
            else {
                this->q2[i][j] = t2[i][j] * (1 - this->tau_s / norm2[i][j]) / (this->mu + this->r);
            }
        }
    } 
}

void AlgorithmChecker::update_lambda() {
    std::cout << std::endl << "============================================================" << std::endl;
    std::cout << "ШАГ 3: Обновление lambda" << std::endl;
    std::cout << "============================================================" << std::endl;

    //double** grad_v1 = allocate_matrix(this->nx - 1, this->ny);
   // double** grad_v2 = allocate_matrix(this->nx, this->ny - 1);
    //gradient(this->v, grad_v1, grad_v2);

    double diff = 0.0;
    double lambda1_new = 0.0,  lambda2_new = 0.0;
    this->diff_max = 0;

    // Обновляем lambda1
    for (int i = 0; i < this->nx - 1; i++) {
        for (int j = 0; j < this->ny; j++) {
            lambda1_new = this->lambda1[i][j] + this->alpha * (grad_v1[i][j] - this->q1[i][j]);
            diff = fabs(lambda1_new - this->lambda1[i][j]);
            if (diff > this->diff_max) this->diff_max = diff;
            this->lambda1[i][j] = lambda1_new;
        }
    }

    // Обновляем lambda2
    for (int i = 0; i < this->nx; i++) {
        for (int j = 0; j < this->ny - 1; j++) {
            lambda2_new = this->lambda2[i][j] + this->alpha * (grad_v2[i][j] - this->q2[i][j]);
            diff = fabs(lambda2_new - this->lambda2[i][j]);
            if (diff > this->diff_max)  this->diff_max = diff;
            this->lambda2[i][j] = lambda2_new;
        }
    }

    calculation_norma(lambda1, lambda2);
    // Сохраняем историю сходимости
    /*
    if (this->convergence_history == nullptr) {
        this->convergence_history = new double[this->current_iteration + 1];
    }
    else {
        double* temp = new double[this->current_iteration + 1];
        memcpy(temp, this->convergence_history, this->current_iteration * sizeof(double));
        delete[] this->convergence_history;
        this->convergence_history = temp;
    }
    this->convergence_history[this->current_iteration - 1] = diff_max;
    */
    std::cout << "Максимальное изменение lambda : " << this->diff_max << std::endl;

   // free_matrix(grad_v1, this->nx - 1);
   // free_matrix(grad_v2, this->nx);
}

void AlgorithmChecker::Save(int bol)
{
    std::cout << "nx = " << nx << " ; " << ny << "\n";
   if (bol == 1)
    { //this->current_iteration % 10 == 0 ||

       
        std::string file_v = this->filename + "output_data/iret" + std::to_string(this->current_iteration) + "_v.txt";
        save_matrix_to_file(this->v, this->nx, this->ny, file_v.c_str());

        /* std::string file_q1 = this->filename + "q1_iteration_" + std::to_string(this->current_iteration) + ".txt";
        save_matrix_to_file(this->q1, this->nx - 1, this->ny, file_q1.c_str());

        std::string file_lambda1 = this->filename + "lambda1_iteration_" + std::to_string(this->current_iteration) + ".txt";
        save_matrix_to_file(this->lambda1, this->nx - 1, this->ny, file_lambda1.c_str());*/

        std::string file_norma1 = this->filename + "output_data/iret" + std::to_string(this->current_iteration) + "_norma11.txt";
        save_matrix_to_file_norma(this->norm1, this->norm2, this->nx - 1, this->ny, file_norma1.c_str());

        /*std::string file_norma1 = this->filename + "output_data/iret" + std::to_string(this->current_iteration) + "_norma12.txt";
        save_matrix_to_file(this->norm2, this->nx - 1, this->ny, file_norma1.c_str());*/

        //std::string file_norm_lambda1 = this->filename + "norm_lambda1_iteration_" + std::to_string(this->current_iteration) + ".txt";
        //save_matrix_to_file(this->norm_lambda1, this->nx - 1, this->ny, file_norm_lambda1.c_str());
    }
   else if (bol == 2)
   {
       std::string file_norma1 = this->filename + "output_data/iret" + std::to_string(this->current_iteration) + "_norma11.txt";
       save_matrix_to_file_norma(this->norm1, this->norm2, this->nx - 1, this->ny, file_norma1.c_str());

     /*  std::string file_norma1 = this->filename + "output_data/iret" + std::to_string(this->current_iteration) + "_norma12.txt";
       save_matrix_to_file(this->norm2, this->nx - 1, this->ny, file_norma1.c_str());*/

        std::ofstream file_diff(this->filename + "diff.txt");
  if (!file_diff.is_open()) {std::cout << "Ошибка открытия файла: " << this->filename + "diff.txt" << std::endl;return;}
  for(int i = 0;i<diff.size();++i)
    file_diff << diff[i] << std::endl;

  file_diff.close();


   }

  /* std::ofstream file_diff(this->filename + "diff.txt", std::ios::app);
   if (!file_diff.is_open()) {std::cout << "Ошибка открытия файла: " << this->filename + "diff.txt" << std::endl;return;}
   file_diff << this->diff_max << std::endl;
   file_diff.close();*/
}
int AlgorithmChecker::run_full_algorithm(int max_iterations, double convergence_tol) {
    std::cout << "ЗАПУСК ПОЛНОГО АЛГОРИТМА" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "Сетка: " << this->nx << " x " << this->ny << " точек" << std::endl;
    std::cout << "Параметры: r=" << this->r << ", mu=" << this->mu << ", tau_s=" << this->tau_s << ", alpha=" << this->alpha << std::endl;
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
        solve_poisson_step(this->f);
        auto poisson_end = std::chrono::high_resolution_clock::now();
        poisson_time = std::chrono::duration<double>(poisson_end - poisson_start).count();

        // Шаг 2: Обновление q с замером времени
        auto update_q_start = std::chrono::high_resolution_clock::now();
        update_q();
        auto update_q_end = std::chrono::high_resolution_clock::now();
        update_q_time = std::chrono::duration<double>(update_q_end - update_q_start).count();

        // Шаг 3: Обновление lambda с замером времени
        auto update_lambda_start = std::chrono::high_resolution_clock::now();
        update_lambda();
        auto update_lambda_end = std::chrono::high_resolution_clock::now();
        update_lambda_time = std::chrono::duration<double>(update_lambda_end - update_lambda_start).count();

       
        // Общее время итерации
        auto iteration_end = std::chrono::high_resolution_clock::now();
        total_time = std::chrono::duration<double>(iteration_end - iteration_start).count();

        // Выводим время выполнения итерации
        print_iteration_time();
        diff.push_back(diff_max);
        if (this->diff_max < convergence_tol) {
            std::cout << std::endl << "СХОДИМОСТЬ ДОСТИГНУТА на итерации " << this->current_iteration << std::endl;
            converged = 1;
            break;
        }

        if (this->current_iteration % 10 == 0) {
            std::cout << std::endl << "--- Прогресс: итерация " << this->current_iteration
                << ", изменение lambda = " << std::scientific << this->diff_max << std::endl;
        }

        if (iteration % 100 == 0)  Save(2);

       // if (iteration + 1 >= max_iterations) Save(1);
        //else Save(3);
        
    }

    Save(1);
    // Вычисляем общее время работы алгоритма
    auto total_end_time = std::chrono::high_resolution_clock::now();
    total_algorithm_time = std::chrono::duration<double>(total_end_time - total_start_time).count();

    // Сохраняем финальные данные
    //save_matrix_to_file(this->v, this->nx, this->ny, "v_matrix_final.txt");
    
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
}//*/




/*
int AlgorithmChecker::run_full_algorithm(int max_iterations, double convergence_tol) {
    std::cout << "ЗАПУСК ПОЛНОГО АЛГОРИТМА" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "Сетка: " << this->nx << " x " << this->ny << " точек" << std::endl;
    std::cout << "Параметры: r =" << this->r << ", mu =" << this->mu << ", tau_s =" << this->tau_s << ", alpha =" << this->alpha << std::endl;
    std::cout << "Максимальное число итераций: " << max_iterations << std::endl;
    std::cout << "Критерий сходимости: " << convergence_tol << std::endl;

    int converged = 0;

    for (int iteration = 0; iteration < max_iterations; iteration++) {
        this->current_iteration = iteration + 1;

        std::cout << std::endl << "============================================================" << std::endl;
        std::cout << "ИТЕРАЦИЯ " << this->current_iteration << std::endl;
        std::cout << "============================================================" << std::endl;

        solve_poisson_step(this->f);
        update_q();
        double diff_max = update_lambda();

        if (diff_max < convergence_tol) {
            std::cout << std::endl << "✓ СХОДИМОСТЬ ДОСТИГНУТА на итерации " << this->current_iteration << std::endl;
            converged = 1;
            break;
        }

        if (this->current_iteration % 10 == 0) {
            std::cout << std::endl << "--- Прогресс: итерация " << this->current_iteration << ", изменение lambda = " << diff_max << std::endl;
        }
    }

    // Сохраняем финальные данные
    save_matrix_to_file(this->v, this->nx, this->ny, "v_matrix_final.txt");

    // ДОБАВИТЬ ФИНАЛЬНУЮ ПРОВЕРКУ СИММЕТРИИ
    std::cout << std::endl << "ФИНАЛЬНАЯ ПРОВЕРКА СИММЕТРИИ:" << std::endl;
    check_symmetry_detailed();

    std::cout << std::endl << "============================================================" << std::endl;
    std::cout << "АЛГОРИТМ ЗАВЕРШЕН" << std::endl;
    if (converged) {
        std::cout << "✓ АЛГОРИТМ СОШЕЛСЯ за " << this->current_iteration << " итераций" << std::endl;
    }
    else {
        std::cout << "✗ ДОСТИГНУТО МАКСИМАЛЬНОЕ ЧИСЛО ИТЕРАЦИЙ (" << max_iterations << ")" << std::endl;
    }

    return converged;
} */