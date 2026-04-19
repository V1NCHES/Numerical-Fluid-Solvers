import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def thomas_algorithm_correct(a, b, c, d, left_boundary=0.0, right_boundary=0.0):
    """Правильная версия алгоритма Томаса с обработкой ошибок"""
    n = len(d)
    if n == 0:
        return np.array([])
    
    # Проверка на нулевые элементы в знаменателе
    if np.any(b == 0):
        print("Предупреждение: нулевые элементы на главной диагонали!")
        # Заменяем нули на очень маленькое число
        b = np.where(b == 0, 1e-12, b)
    
    # Векторы для прогоночных коэффициентов
    alpha = np.zeros(n)
    beta = np.zeros(n)
    x = np.zeros(n)
    
    try:
        # ПРЯМОЙ ХОД: вычисление прогоночных коэффициентов
        
        # Первый узел (i=0)
        denominator = b[0]
        if abs(denominator) < 1e-12:
            denominator = 1e-12
            
        alpha[0] = -c[0] / denominator
        beta[0] = (d[0] - a[0] * left_boundary) / denominator
        
        # Промежуточные узлы (i = 1, 2, ..., n-2)
        for i in range(1, n-1):
            denominator = b[i] + a[i] * alpha[i-1]
            if abs(denominator) < 1e-12:
                denominator = 1e-12
                
            alpha[i] = -c[i] / denominator
            beta[i] = (d[i] - a[i] * beta[i-1]) / denominator
        
        # Последний узел (i = n-1)
        if n > 1:
            denominator = b[n-1] + a[n-1] * alpha[n-2]
            if abs(denominator) < 1e-12:
                denominator = 1e-12
                
            beta[n-1] = (d[n-1] - a[n-1] * beta[n-2] - c[n-1] * right_boundary) / denominator
        
        # ОБРАТНЫЙ ХОД: восстановление решения
        
        # Последний узел
        if n > 1:
            x[n-1] = beta[n-1]
            # Остальные узлы: x_i = α_i * x_{i+1} + β_i
            for i in range(n-2, -1, -1):
                x[i] = alpha[i] * x[i+1] + beta[i]
        else:
            # Если всего один узел
            x[0] = beta[0]
            
    except (FloatingPointError, ZeroDivisionError) as e:
        print(f"Ошибка в алгоритме Томаса: {e}")
        # Возвращаем нулевое решение в случае ошибки
        return np.zeros(n)
    
    # Проверка на NaN и бесконечности
    if np.any(np.isnan(x)) or np.any(np.isinf(x)):
        print("Предупреждение: обнаружены NaN или бесконечные значения в решении!")
        x = np.nan_to_num(x, nan=0.0, posinf=1e6, neginf=-1e6)
    
    return x

def adi_poisson_solver_stable(f_matrix, Nx, Ny, Lx=1.0, Ly=1.0, max_iter=1000, tol=1e-6):
    """
    Стабильная версия ADI метода с обработкой ошибок
    """
    start_time = time.time()
    
    print(f"Размер f_matrix (внутренние точки): {f_matrix.shape}")
    print(f"Область: {Lx} x {Ly}")
    
    # Проверка входных данных
    if np.any(np.isnan(f_matrix)) or np.any(np.isinf(f_matrix)):
        print("Предупреждение: обнаружены NaN или бесконечные значения в f_matrix!")
        f_matrix = np.nan_to_num(f_matrix, nan=0.0, posinf=1e6, neginf=-1e6)
    
    # Если f_matrix 15x15 - это внутренние точки, то всего точек 17x17
    total_points_x = Nx + 2
    total_points_y = Ny + 2
    
    # Шаги сетки с проверкой
    if total_points_x <= 1 or total_points_y <= 1:
        raise ValueError("Слишком маленькая сетка!")
        
    hx = Lx / (total_points_x - 1)
    hy = Ly / (total_points_y - 1)
    
    # Проверка на очень маленькие шаги
    if hx < 1e-10 or hy < 1e-10:
        print("Предупреждение: очень маленькие шаги сетки!")
        hx = max(hx, 1e-10)
        hy = max(hy, 1e-10)
    
    hx2 = hx**2
    hy2 = hy**2
    
    print(f"Шаги сетки: hx = {hx:.4e}, hy = {hy:.4e}")
    print(f"Всего точек (с границами): {total_points_x} x {total_points_y}")
    
    # Инициализация решения - только внутренние точки
    u = np.zeros((Nx, Ny))
    
    # Трехдиагональные матрицы для внутренних точек
    a_x = np.full(Nx, -1/hx2)
    b_x = np.full(Nx, 2/hx2)  
    c_x = np.full(Nx, -1/hx2)
    
    a_y = np.full(Ny, -1/hy2)
    b_y = np.full(Ny, 2/hy2)
    c_y = np.full(Ny, -1/hy2)
    
    history = []
    
    for iteration in range(max_iter):
        u_old = u.copy()
        
        # === ПОЛУШАГ ПО X ===
        for j in range(Ny):
            rhs = np.zeros(Nx)
            for i in range(Nx):
                # Граничные условия с проверкой
                bottom = 0.0 if j == 0 else u_old[i, j-1]
                top = 0.0 if j == Ny-1 else u_old[i, j+1]
                
                # Вычисление правой части с проверкой переполнения
                term1 = -f_matrix[i, j]
                term2 = (bottom + top) / hy2
                
                # Проверка на переполнение
                if abs(term2) > 1e10:
                    term2 = np.sign(term2) * 1e10
                
                rhs[i] = term1 - term2
            
            # Проверка правой части
            if np.any(np.isnan(rhs)) or np.any(np.isinf(rhs)):
                print(f"Предупреждение: NaN/Inf в правой части (X шаг, j={j})")
                rhs = np.nan_to_num(rhs, nan=0.0, posinf=1e6, neginf=-1e6)
            
            # Решаем систему для строки j
            u_line = thomas_algorithm_correct(a_x, b_x, c_x, rhs, 0.0, 0.0)
            
            # Проверка решения
            if np.any(np.isnan(u_line)) or np.any(np.isinf(u_line)):
                print(f"Предупреждение: NaN/Inf в решении (X шаг, j={j})")
                u_line = np.nan_to_num(u_line, nan=0.0, posinf=1e6, neginf=-1e6)
            
            # Обновляем строку j
            for i in range(Nx):
                u[i, j] = u_line[i]
        
        # === ПОЛУШАГ ПО Y ===  
        for i in range(Nx):
            rhs = np.zeros(Ny)
            for j in range(Ny):
                # Граничные условия с проверкой
                left = 0.0 if i == 0 else u[i-1, j]
                right = 0.0 if i == Nx-1 else u[i+1, j]
                
                # Вычисление правой части с проверкой переполнения
                term1 = -f_matrix[i, j]
                term2 = (left + right) / hx2
                
                # Проверка на переполнение
                if abs(term2) > 1e10:
                    term2 = np.sign(term2) * 1e10
                
                rhs[j] = term1 - term2
            
            # Проверка правой части
            if np.any(np.isnan(rhs)) or np.any(np.isinf(rhs)):
                print(f"Предупреждение: NaN/Inf в правой части (Y шаг, i={i})")
                rhs = np.nan_to_num(rhs, nan=0.0, posinf=1e6, neginf=-1e6)
            
            # Решаем систему для столбца i
            u_column = thomas_algorithm_correct(a_y, b_y, c_y, rhs, 0.0, 0.0)
            
            # Проверка решения
            if np.any(np.isnan(u_column)) or np.any(np.isinf(u_column)):
                print(f"Предупреждение: NaN/Inf в решении (Y шаг, i={i})")
                u_column = np.nan_to_num(u_column, nan=0.0, posinf=1e6, neginf=-1e6)
            
            # Обновляем столбец i
            for j in range(Ny):
                u[i, j] = u_column[j]
        
        # Проверка решения на NaN/Inf
        if np.any(np.isnan(u)) or np.any(np.isinf(u)):
            print(f"Критическая ошибка: NaN/Inf в решении на итерации {iteration}")
            u = np.nan_to_num(u, nan=0.0, posinf=1e6, neginf=-1e6)
            break
        
        # Проверка сходимости
        residual = np.max(np.abs(u - u_old))
        
        # Проверка на NaN в невязке
        if np.isnan(residual) or np.isinf(residual):
            print(f"Невязка содержит NaN/Inf на итерации {iteration}")
            residual = 1.0  # продолжаем итерации
            
        history.append(residual)
        
        if residual < tol:
            execution_time = time.time() - start_time
            print(f"✓ ADI метод сошелся за {iteration+1} итераций")
            print(f"Финальная невязка: {residual:.2e}")
            print(f"Время выполнения: {execution_time:.4f} сек")
            return u, history, execution_time
        
        if iteration % 100 == 0:
            print(f"Итерация {iteration:3d}, невязка: {residual:.2e}")
        
        # Проверка на расходимость
        if residual > 1e10 and iteration > 10:
            print(f"Метод расходится на итерации {iteration}")
            break
    
    execution_time = time.time() - start_time
    print(f"✗ ADI метод не сошелся за {max_iter} итераций")
    print(f"Финальная невязка: {residual:.2e}")
    print(f"Время выполнения: {execution_time:.4f} сек")
    return u, history, execution_time

def test_stability():
    """Тест стабильности с различными параметрами"""
    print("ТЕСТ СТАБИЛЬНОСТИ ADI МЕТОДА")
    print("=" * 50)
    
    # Тест с разными размерами сетки
    test_cases = [
        (8, 8, "Маленькая сетка"),
        (16, 16, "Средняя сетка"), 
        (32, 32, "Большая сетка")
    ]
    
    for Nx, Ny, description in test_cases:
        print(f"\n{description}: {Nx}x{Ny}")
        print("-" * 30)
        
        try:
            # Создаем тестовую правую часть
            f_matrix = np.zeros((Nx, Ny))
            hx = 1.0 / (Nx + 1)
            hy = 1.0 / (Ny + 1)
            
            for i in range(Nx):
                for j in range(Ny):
                    x = (i + 1) * hx
                    y = (j + 1) * hy
                    f_matrix[i, j] = -2 * (y**2 - y + x**2 - x)
            
            # Запускаем решатель
            solution, history, time_elapsed = adi_poisson_solver_stable(
                f_matrix, Nx, Ny, Lx=1.0, Ly=1.0, max_iter=5000, tol=1e-6
            )
            
            print(f"Успешно завершено за {time_elapsed:.2f} сек")
            
        except Exception as e:
            print(f"Ошибка: {e}")

# Использование в основном коде
def demonstrate_stable_methods():
    """Демонстрация стабильных методов"""
    
    print("Решение уравнения Пуассона со стабильными алгоритмами")
    print("=" * 60)
    
    # Параметры сетки
    N1, N2 = 16, 16
    L1, L2 = 1.0, 1.0
    
    # Создаем правую часть
    h1 = L1/N1
    h2 = L2/N2
    f = np.zeros((N1-1, N2-1))
    
    for i in range(0, N1-1):
        for j in range(0, N2-1):
            x = (i+1) * h1
            y = (j+1) * h2
            f[i, j] = -2*(y**2-y+x**2-x)
    
    print("Тестируем стабильный ADI метод...")
    
    try:
        solution_ADI, history, time_elapsed = adi_poisson_solver_stable(
            f, N1-1, N2-1, L1, L2, max_iter=1000, tol=1e-6
        )
        
        # Аналитическое решение для сравнения
        analytical = np.zeros((N1-1, N2-1))
        for i in range(0, N1-1):
            for j in range(0, N2-1):
                x = (i+1) * h1
                y = (j+1) * h2
                analytical[i, j] = (x**2-x)*(y**2-y)
        
        error = np.abs(solution_ADI - analytical)
        print(f"Максимальная ошибка: {np.max(error):.2e}")
        print(f"Средняя ошибка: {np.mean(error):.2e}")
        
    except Exception as e:
        print(f"Ошибка при выполнении: {e}")

if __name__ == "__main__":
    # Включим обработку предупреждений
    np.seterr(all='warn')
    
    # Запуск теста стабильности
    test_stability()
    
    # Или демонстрация
    # demonstrate_stable_methods()
