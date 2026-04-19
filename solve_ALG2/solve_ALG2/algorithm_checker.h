#pragma once
#include <chrono> 
#include <vector>
#include <iostream>
#include <locale>
#include <windows.h>
//class CorrectPoissonSolver2D;
class PoissonSolver2D;


class AlgorithmChecker {
public:
    // Параметры сетки и области
    int nx;          // Количество узлов по оси X
    int ny;          // Количество узлов по оси Y
    double Lx;       // Длина области по оси X
    double Ly;       // Длина области по оси Y
    double hx;       // Шаг сетки по оси X
    double hy;       // Шаг сетки по оси Y

    // Матрицы для алгоритма (двойные указатели для динамических 2D массивов)
    double** lambda1;  // Множители Лагранжа для компоненты X (размер: (nx-1) x ny)
    double** lambda2;  // Множители Лагранжа для компоненты Y (размер: nx x (ny-1))
    double** q1;       // Вспомогательные переменные для компоненты X (размер: (nx-1) x ny)
    double** q2;       // Вспомогательные переменные для компоненты Y (размер: nx x (ny-1))
    double** v;        // Основная переменная (потенциал) (размер: nx x ny)
    double** v_f;      // Основная переменная для вычесления во внутреней сеткке скорости (размер: (nx -2) x (ny - 2) )
    double** f_i;       // Правая часть во внутреней сеткке (размер: (nx -2) x (ny - 2) )
    double** div_lambda; // Дивергенция во внутреней сеткке Ламбды (размер: (nx -2) x (ny - 2) )
    double** div_q;      // Дивергенция во внутреней сеткке q (размер: (nx -2) x (ny - 2) )
    double** F_i; // Матрица правых частей с граничными условиями nx*ny

     // Матрицы для вычесления новых q
    double** grad_v1;//
    double** grad_v2;//
    double** t1;
    double** t2;
    double** t1_extended;
    double** t2_extended;
    double** t1_avg;
    double** t2_avg;
    double** norm1;
    double** norm2;
    double** norm_lambda1;

    // Параметры алгоритма
    double r;           // Параметр регуляризации
    double tau_s;       // Параметр сдвига (пороговое значение)
    double mu;          // Параметр вязкости/трения
    double alpha;       // Параметр скорости обучения/релаксации
    int current_iteration;  // Текущая итерация алгоритма
    double diff_max; // сходимость алгоритма

    // Правая часть уравнения и решатель
    double** f;                    // Перепад давления (размер: (nx-2) x (ny-2))
    //CorrectPoissonSolver2D* poisson_solver;  // Указатель на решатель уравнения Пуассона
    PoissonSolver2D* poisson_solver;  // Указатель на решатель уравнения Пуассона

    // Координатные оси
    double* x;          // Массив координат по оси X
    double* y;          // Массив координат по оси Y

    // История сходимости и состояний
    double* convergence_history;  // Массив для хранения истории изменений (сходимости)
    double*** v_history;          // Трехмерный массив для хранения истории состояний v
    int history_size;             // Текущий размер истории

    //double*** q1_history;

     // Добавляем поля для измерения времени
    double total_time;
    double poisson_time;
    double update_q_time;
    double update_lambda_time;

    //файл для записи данные
    std::string filename;

    //Вектора для вычесления скорость ядра v_plug /  DLT - это показатель сходимости численного решения
    std::vector <double> v_plug;
    std::vector <double> DLT;
    std::vector <double> diff;
    // Методы для измерения времени
    void reset_timers();
    void print_iteration_time();

    // Конструктор: инициализирует все параметры и выделяет память
    AlgorithmChecker(double hx, double hy, double Lx, double Ly, double r, double mu, double tau_s, double alpha, const std::string& filename);

    //AlgorithmChecker(int nx = 9, int ny = 9, double r = 1.0, double mu = 1.0, double tau_s = 0.2, double alpha = 1.0, const std::string & filename = "data/iterator_");
    

    // Деструктор: освобождает всю выделенную память
    ~AlgorithmChecker();

    // Методы управления данными
    void create_f_matrix();  // Создает тестовую матрицу правой части f
    void save_matrix_to_file(double** matrix, int rows, int cols, const char* filename);  // Сохраняет матрицу в файл
    static double** allocate_matrix(int rows, int cols);  // Выделяет память для 2D матрицы
    static void free_matrix(double** matrix, int rows);   // Освобождает память 2D матрицы

    // Математические операции
    void divergence(double** field1, double** field2, double** div);  // Вычисляет дивергенцию поля
    void gradient(double** v, double** grad_v1, double** grad_v2);    // Вычисляет градиент скалярного поля

    // Основные шаги алгоритма
    int solve_poisson_step(double** f);  // Шаг 1: Решает уравнение Пуассона для v
    void update_q();                     // Шаг 2: Обновляет вспомогательные переменные q
    void update_lambda();              // Шаг 3: Обновляет множители Лагранжа lambda

    // Главный метод алгоритма
    int run_full_algorithm(int max_iterations = 1000, double convergence_tol = 1e-6);  // Запускает полный алгоритм

    // Метод проверки и диагностики
    void check_symmetry_detailed();  // Проверяет симметричность матрицы v для диагностики

    void Save(int bol); //Сохраняем данны на кажом шаге или 

    void calculation_norma(double** Q1, double** Q2);
    void save_matrix_to_file_norma(double** norna1, double** norma2, int rows, int cols, const char* filename);
};//*/


void main_Algoritm();
