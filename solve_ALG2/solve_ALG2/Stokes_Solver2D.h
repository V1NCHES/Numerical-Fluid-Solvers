#pragma once
#include <chrono> 
#include <vector>
#include <iostream>



class StokesSolver2D {
public:
    // Параметры сетки и области
    int nx1, nx2;          // Количество клеток , всего ячеек nx1 + 1 + nx2
    int ny1, ny2;          // Количество клеток ,  всего ячеек ny1 + 1 + ny2
    double nx_d;           // Количесво ячеек для дырки
    double Lx1, Lx2;       // Длина области по оси X
    double Ly1, Ly2;       // Длина области по оси Y

    double hx;       // Шаг сетки по оси X
    double hy;       // Шаг сетки по оси Y

    double d_t, h2, sm2, osp2; //  ПТМ константы

    std::vector<std::vector<double>> v1;        // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2) x ( ny1 + 1 )
    std::vector<std::vector<double>> v2;        // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)
    std::vector<std::vector<double>> u1;        // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2) x ( ny1 + 1 )
    std::vector<std::vector<double>> u2;        // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)

    std::vector<std::vector<double>> fv1;        // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2) x ( ny1 + 1 )
    std::vector<std::vector<double>> fv2;        // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)
    std::vector<std::vector<double>> fu1;        // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2) x ( ny1 + 1 )
    std::vector<std::vector<double>> fu2;        // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)


    std::vector<std::vector<double>> Fv1;        // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2) x ( ny1 + 1 )
    std::vector<std::vector<double>> Fv2;        // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)
    std::vector<std::vector<double>> Fu1;        // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2) x ( ny1 + 1 )
    std::vector<std::vector<double>> Fu2;        // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)

    std::vector<std::vector<double>> p1;        // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2) x ( ny1 + 1 )
    std::vector<std::vector<double>> p2;        // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)

    std::vector<std::vector<double>> g1;        // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2) x ( ny1 + 1 )
    std::vector<std::vector<double>> g2;        // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)

    std::vector<std::vector<double>> r1;        // Основная переменная (потенциал) (размер: (nx1 + 1 +  nx2) x ( ny1 + 1 )
    std::vector<std::vector<double>> r2;        // Основная переменная (потенциал) (размер: (nx1 + 1) x (ny2)

    //// Матрицы для алгоритма (двойные указатели для динамических 2D массивов)
    //std::vector<std::vector<double>> lambda11;
    //std::vector<std::vector<double>> lambda12;
    //std::vector<std::vector<double>> lambda21;
    //std::vector<std::vector<double>> lambda22;

    //std::vector<std::vector<double>> q11;
    //std::vector<std::vector<double>> q12;
    //std::vector<std::vector<double>> q21;
    //std::vector<std::vector<double>> q22;

    //std::vector<std::vector<double>> div_lambda1;
    //std::vector<std::vector<double>> div_lambda2;

    //std::vector<std::vector<double>> div_q1;
    //std::vector<std::vector<double>> div_q2;


    //// Матрицы для вычесления новых q
    //std::vector<std::vector<double>> grad_v11;
    //std::vector<std::vector<double>> grad_v12;
    //std::vector<std::vector<double>> grad_v21;
    //std::vector<std::vector<double>> grad_v22;

    //std::vector<std::vector<double>> t11;
    //std::vector<std::vector<double>> t12;
    //std::vector<std::vector<double>> t21;
    //std::vector<std::vector<double>> t22;

    //std::vector<std::vector<double>> t11_extended;
    //std::vector<std::vector<double>> t12_extended;
    //std::vector<std::vector<double>> t21_extended;
    //std::vector<std::vector<double>> t22_extended;

    //std::vector<std::vector<double>> t11_avg;
    //std::vector<std::vector<double>> t12_avg;
    //std::vector<std::vector<double>> t21_avg;
    //std::vector<std::vector<double>> t22_avg;

    //std::vector<std::vector<double>> norma11;
    //std::vector<std::vector<double>> norma12;
    //std::vector<std::vector<double>> norma21;
    //std::vector<std::vector<double>> norma22;

    //std::vector<std::vector<double>> norma_lambda11;
    //std::vector<std::vector<double>> norma_lambda21;

    // Параметры алгоритма
    double r;           // Параметр регуляризации
    double tau_s;       // Параметр сдвига (пороговое значение)
    double mu;          // Параметр вязкости/трения
    double alpha;       // Параметр скорости обучения/релаксации
    int current_iteration;  // Текущая итерация алгоритма
    double diff_max; // сходимость алгоритма

     // Добавляем поля для измерения времени
    double total_time;
    double poisson_time;
    double update_q_time;
    double update_lambda_time;

    //файл для записи данные
    std::string filename;

    // Методы для измерения времени
    void reset_timers() { total_time = 0.0;  poisson_time = 0.0; update_q_time = 0.0; update_lambda_time = 0.0; }
    void print_iteration_time();

    // Конструктор: инициализирует все параметры и выделяет память
    StokesSolver2D(double hx, double hy, double Lx1, double Lx2, double Ly1, double Ly2, double r, double mu, double tau_s, double alpha, const std::string& filename);
    StokesSolver2D();
    // Деструктор: освобождает всю выделенную память
    ~StokesSolver2D();

    // Главный метод алгоритма
    int run_full_algorithm(int max_iterations = 1000, double convergence_tol = 1e-6);  // Запускает полный алгоритм


    // Методы управления данными
    void create_initial_conditions();  // Задаем начальные усвлоия для v, u, fv, fu

    // Математические операции
    //void divergence(std::vector<std::vector<double>>& field11, std::vector<std::vector<double>>& field12, std::vector<std::vector<double>>& div1, std::vector<std::vector<double>>& field21, std::vector<std::vector<double>>& field22, std::vector<std::vector<double>>& div2);  // Вычисляет дивергенцию поля
    //void gradient_v();

    // Основные шаги алгоритма
    void update_Poisson(); 
    void Residual_calculation(); // g = Bx u + By v
    void Two_step_iterative_method(); // p *  <- Ap = g, A = Bx + By

    //Шаг 0: Обновлаем уранение Пуассона F(x,y) = (f + div_lambda - r * div_q) / r; 
   // void update_Stopkes(std::vector<std::vector<double>>& F1, std::vector<std::vector<double>>& F2, std::vector<std::vector<double>>& f1, std::vector<std::vector<double>>& f2);
    //Шаг 0.1: Обновлаем Стоксовый член  F(x,y) -= BxP; (ByP в зависемости какой v или ux`
    // Решаем ур Пуассона МПТ  
    void PoissonSolver2D(double EPS, int max_iter, std::vector<std::vector<double>>& V1, std::vector<std::vector<double>>& V2, std::vector<std::vector<double>>& F1, std::vector<std::vector<double>>& F2);
   //Шаг 1: Решает уравнение Пуассона для v
    double U_ij(double u_ij, double u_ip1j, double u_ijp1, double u_im1j, double u_ijm1, double f_ij);
    double Error(double u_ij, double u_ip1j, double u_ijp1, double u_im1j, double u_ijm1, double f_ij);
    double Error_Area(const std::vector<std::vector<double>>& data1, const std::vector<std::vector<double>>& data2);

    //void update_q();                     // Шаг 2: Обновляет вспомогательные переменные q
    //void calculation_t(std::vector<std::vector<double>>& lamda, std::vector<std::vector<double>>& grad_v, std::vector<std::vector<double>>& t, std::vector<std::vector<double>>& t_extended, int Nx, int Ny, int k, int l);
    //void averaging_t(std::vector<std::vector<double>>& t_avg, std::vector<std::vector<double>>& t_extended, int Nx, int Ny);
    //void calculation_norma(std::vector<std::vector<double>>& t, std::vector<std::vector<double>>& t_avg, std::vector<std::vector<double>>& norma, int Nx, int Ny);
    //void calculation_q(std::vector<std::vector<double>>& t, std::vector<std::vector<double>>& norma, std::vector<std::vector<double>>& q, int Nx, int Ny);

    //void update_lambda();              // Шаг 3: Обновляет множители Лагранжа lambda
    //void calculation_lamda(std::vector<std::vector<double>>& lambda, std::vector<std::vector<double>>& grad_v, std::vector<std::vector<double>>& q, int Nx, int Ny);





    void Save(); //Сохраняем данны на кажом шаге или 
    void save_V_file(const std::vector<std::vector<double>>& data1, const std::vector<std::vector<double>>& data2, int Nx1, int Ny1, int Nx2, int Ny2, const std::string& filename);
    void check_symmetry_detailed();  // Проверяет симметричность матрицы v для диагностики

};//*/

//void print_v(std::vector<std::vector<double>>& print1, std::vector<std::vector<double>>& print2, int nx1, int nx2, int ny1, int ny2);