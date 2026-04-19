import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def visualize_solutions(numerical, analytical, error, L1, L2):
    """Визуализация численного и аналитического решений"""
    
    N1_total, N2_total = numerical.shape
    x = np.linspace(0, L1, N1_total)
    y = np.linspace(0, L2, N2_total)
    X, Y = np.meshgrid(x, y, indexing='ij')
    '''
    fig = plt.figure(figsize=(18, 5))
    
    # Численное решение
    ax1 = fig.add_subplot(131, projection='3d')
    surf1 = ax1.plot_surface(X, Y, numerical, cmap='viridis', alpha=0.8)
    ax1.set_title('Численное решение\n(по точным формулам файла)')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    ax1.set_zlabel('u(x,y)')
    
    # Аналитическое решение
    ax2 = fig.add_subplot(132, projection='3d')
    surf2 = ax2.plot_surface(X, Y, analytical, cmap='viridis', alpha=0.8)
    ax2.set_title('Аналитическое решение')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_zlabel('u(x,y)')
    # Погрешность
    ax3 = fig.add_subplot(133, projection='3d')
    surf3 = ax3.plot_surface(X, Y, error, cmap='hot', alpha=0.8)
    ax3.set_title('Погрешность')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('|error|')
    
    plt.tight_layout()
    plt.show()
    '''
    # Дополнительно: контурные графики
    fig2, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    im1 = axes[0].contourf(X, Y, numerical, levels=20, cmap='viridis')
    axes[0].set_title('Численное решение')
    axes[0].set_xlabel('X')
    axes[0].set_ylabel('Y')
    plt.colorbar(im1, ax=axes[0])
    
    im2 = axes[1].contourf(X, Y, analytical, levels=20, cmap='viridis')
    axes[1].set_title('Аналитическое решение')
    axes[1].set_xlabel('X')
    axes[1].set_ylabel('Y')
    plt.colorbar(im2, ax=axes[1])
    
    im3 = axes[2].contourf(X, Y, error, levels=20, cmap='hot')
    axes[2].set_title('Погрешность')
    axes[2].set_xlabel('X')
    axes[2].set_ylabel('Y')
    plt.colorbar(im3, ax=axes[2])
    
    plt.tight_layout()
    plt.show()
# =============================================================================
# Преобразование фурье DST
# =============================================================================

def lamda_k_0(k,h,N):   
    return 4*(np.sin(k*np.pi/(2*N))**2)/(h**2)


def phi_k2(f_ij,i,k2,N2):
    SUM = 0.0
    for j in range(1,N2,1):
        SUM += f_ij[i-1][j-1]*np.sin(k2*np.pi*j/N2)
        
    return SUM

def phi_k1_k2(f_ij,k1,N1,k2,N2):   
    SUM = 0.0
    for i in range(1,N1,1):
        SUM += phi_k2(f_ij, i,k2,N2)*np.sin(k1*np.pi*i/N1)
        
    return SUM

def u_k2(f_ij, i, h1,N1,k2,h2,N2):   
    SUM = 0.0
    for k1 in range(1,N1,1):
        SUM += phi_k1_k2(f_ij,k1,N1,k2,N2)*np.sin(k1*np.pi*i/N1)/(lamda_k_0(k1,h1,N1)+ lamda_k_0(k2,h2,N2) )       
    return SUM

def u_ij(f_ij,i,j, h1,h2, N1,N2):
    
    SUM = 0.0
    for k2 in range(1,N1,1):
        SUM += u_k2(f_ij, i, h1,N1,k2,h2,N2)*np.sin(k2*np.pi*j/N2)        
    return SUM

def lamda_k_0_correct(k, h, L):
    return 4/(h**2) * np.sin(np.pi * k * h / (2 * L))**2

def solve_u_DST(f_inner, N1, N2, h1, h2, L1, L2):
    start_time = time.time()
    u = np.zeros((N1-1, N2-1))
    print(N1-1,N2-1)
    for i in range(1, N1):
        for j in range(1, N2):
            u[i-1,j-1] = u_ij(f_inner,i,j, h1,h2, N1,N2)*4/(N1*N2)
    execution_time = time.time() - start_time
    print(f"Финальное Время выполнения DST: {execution_time:.4f} сек")    
    return u
    
# =============================================================================
# Метод переменных направлений (ADI)
# =============================================================================

def thomas_algorithm_correct(a, b, c, d, left_boundary=0.0, right_boundary=0.0):
    """Правильная версия алгоритма Томаса с прогоночными коэффициентами"""
    n = len(d)
    if n == 0:
        return np.array([])
    
    # Создаем копии для безопасности
    a = a.copy()
    b = b.copy()
    c = c.copy()
    d = d.copy()
    
    # Векторы для прогоночных коэффициентов
    alpha = np.zeros(n)
    beta = np.zeros(n)
    x = np.zeros(n)
    
    # ПРЯМОЙ ХОД: вычисление прогоночных коэффициентов
    
    # Первый узел (i=0)
    # α₀ = -c₀ / b₀
    # β₀ = (d₀ - a₀ * u_left) / b₀
    alpha[0] = -c[0] / b[0]
    beta[0] = (d[0] - a[0] * left_boundary) / b[0]
    
    # Промежуточные узлы (i = 1, 2, ..., n-2)
    for i in range(1, n-1):
        denominator = b[i] + a[i] * alpha[i-1]
        alpha[i] = -c[i] / denominator
        beta[i] = (d[i] - a[i] * beta[i-1]) / denominator
    
    # Последний узел (i = n-1)
    if n > 1:
        denominator = b[n-1] + a[n-1] * alpha[n-2]
        # Учитываем правую границу в последнем уравнении
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
    
    return x

def adi_poisson_solver(f_matrix, Nx, Ny, Lx=1.0, Ly=1.0, max_iter=1000, tol=1e-6):
    """
    ADI метод с правильным алгоритмом Томаса
    f_matrix - матрица БЕЗ граничных точек
    """
    start_time = time.time()
    
    print(f"Размер f_matrix (внутренние точки): {f_matrix.shape}")
    print(f"Область: {Lx} x {Ly}")
    
    # Если f_matrix 15x15 - это внутренние точки, то всего точек 17x17
    # (15 внутренних + 2 границы)
    total_points_x = Nx + 2
    total_points_y = Ny + 2
    
    # Шаги сетки
    hx = Lx / (total_points_x - 1)
    hy = Ly / (total_points_y - 1)
    hx2 = hx**2
    hy2 = hy**2
    
    print(f"Шаги сетки: hx = {hx:.4f}, hy = {hy:.4f}")
    print(f"Всего точек (с границами): {total_points_x} x {total_points_y}")
    
    # Инициализация решения - только внутренние точки
    u = np.zeros((Nx, Ny))  # 15x15 - только внутренние
    
    # Трехдиагональные матрицы для внутренних точек
    # Коэффициенты для системы: a[i]*u[i-1] + b[i]*u[i] + c[i]*u[i+1] = d[i]
    a_x = np.full(Nx, -1/hx2)  # коэффициент перед u[i-1,j]
    b_x = np.full(Nx, 2/hx2)   # коэффициент перед u[i,j]  
    c_x = np.full(Nx, -1/hx2)  # коэффициент перед u[i+1,j]
    
    a_y = np.full(Ny, -1/hy2)  # коэффициент перед u[i,j-1]
    b_y = np.full(Ny, 2/hy2)   # коэффициент перед u[i,j]
    c_y = np.full(Ny, -1/hy2)  # коэффициент перед u[i,j+1]
    
    history = []
    
    for iteration in range(max_iter):
        u_old = u.copy()
        
        # === ПОЛУШАГ ПО X ===
        # Решаем системы вдоль строк (фиксируем j, решаем по i)
        for j in range(Ny):  # j от 0 до Ny-1 (внутренние точки по y)
            rhs = np.zeros(Nx)
            for i in range(Nx):  # i от 0 до Nx-1 (внутренние точки по x)
                # Граничные условия: если точка на краю внутренней области,
                # то соседняя точка - это граница (u = 0)
                bottom = 0.0 if j == 0 else u_old[i, j-1]      # нижняя граница
                top = 0.0 if j == Ny-1 else u_old[i, j+1]      # верхняя граница
                
                # Правая часть уравнения для полушага по X:
                # -f - (u_верх + u_низ) / hy²
                rhs[i] = -f_matrix[i, j] - (bottom + top) / hy2
            
            # Граничные условия для прогонки по X
            # На левой границе (i = -1): u = 0
            # На правой границе (i = Nx): u = 0  
            left_boundary_x = 0.0
            right_boundary_x = 0.0
            
            # Решаем трехдиагональную систему для строки j
            u_line = thomas_algorithm_correct(a_x, b_x, c_x, rhs, left_boundary_x, right_boundary_x)
            
            # Обновляем строку j промежуточного решения u*
            for i in range(Nx):
                u[i, j] = u_line[i]
        
        # === ПОЛУШАГ ПО Y ===
        # Решаем системы вдоль столбцов (фиксируем i, решаем по j)  
        for i in range(Nx):  # i от 0 до Nx-1 (внутренние точки по x)
            rhs = np.zeros(Ny)
            for j in range(Ny):  # j от 0 до Ny-1 (внутренние точки по y)
                # Граничные условия
                left = 0.0 if i == 0 else u[i-1, j]       # левая граница
                right = 0.0 if i == Nx-1 else u[i+1, j]   # правая граница
                
                # Правая часть уравнения для полушага по Y:
                # -f - (u_лево + u_право) / hx²
                rhs[j] = -f_matrix[i, j] - (left + right) / hx2
            
            # Граничные условия для прогонки по Y
            # На нижней границе (j = -1): u = 0
            # На верхней границе (j = Ny): u = 0
            bottom_boundary_y = 0.0
            top_boundary_y = 0.0
            
            # Решаем трехдиагональную систему для столбца i
            u_column = thomas_algorithm_correct(a_y, b_y, c_y, rhs, bottom_boundary_y, top_boundary_y)
            
            # Обновляем столбец i окончательного решения u^{n+1}
            for j in range(Ny):
                u[i, j] = u_column[j]
        
        # Проверка сходимости
        residual = np.max(np.abs(u - u_old))
        history.append(residual)
        
        if residual < tol:
            execution_time = time.time() - start_time
            print(f"✓ ADI метод сошелся за {iteration+1} итераций")
            print(f"Финальная невязка: {residual:.2e}")
            print(f"Время выполнения: {execution_time:.4f} сек")
            return u#, history, execution_time
        
        if iteration % 100 == 0:
            print(f"Итерация {iteration:3d}, невязка: {residual:.2e}")
    
    execution_time = time.time() - start_time
    print(f"✗ ADI метод не сошелся за {max_iter} итераций")
    print(f"Финальная невязка: {residual:.2e}")
    print(f"Время выполнения: {execution_time:.4f} сек")
    return u#, history, execution_time

# =============================================================================
#     Попеременно-треугольный метод для уравнения Пуассона
# =============================================================================
def alternating_triangle_method_poisson_solver(f_matrix, Nx, Ny, Lx=1.0, Ly=1.0, max_iter=1000, tol=1e-6):
    """
    Попеременно-треугольный метод для уравнения Пуассона
    """
    start_time = time.time()
    
    # Раздельные шаги сетки
    total_points_x = Nx + 2
    total_points_y = Ny + 2
    hx = Lx / (total_points_x - 1)
    hy = Ly / (total_points_y - 1)
    hx2 = hx**2
    hy2 = hy**2
    
    # Параметр релаксации для прямоугольной области
    tau = 0.8 / (2/hx2 + 2/hy2)  # оптимальный для прямоугольника
    
    print(f"Улучшенный попеременно-треугольный метод:")
    print(f"Область: {Lx} x {Ly}")
    print(f"Внутренние точки: {Nx} x {Ny}")
    print(f"Шаги: hx = {hx:.4f}, hy = {hy:.4f}")
    print(f"Параметр релаксации: τ = {tau:.2e}")
    
    u = np.zeros((Nx, Ny))
    history = []
    
    for iteration in range(max_iter):
        u_old = u.copy()
        
        # Прямой ход
        for i in range(Nx):
            for j in range(Ny):
                if (i + j) % 2 == 0:
                    # Более точная аппроксимация для прямоугольника
                    left = 0.0 if i == 0 else u[i-1, j]
                    right = 0.0 if i == Nx-1 else u[i+1, j]
                    bottom = 0.0 if j == 0 else u[i, j-1]
                    top = 0.0 if j == Ny-1 else u[i, j+1]
                    
                    # Оператор Лапласа с раздельными шагами
                    laplacian_u = (left + right - 2*u[i, j]) / hx2 + (bottom + top - 2*u[i, j]) / hy2
                    
                    u[i, j] += tau * (f_matrix[i, j] + laplacian_u)
        
        # Обратный ход
        for i in range(Nx-1, -1, -1):
            for j in range(Ny-1, -1, -1):
                if (i + j) % 2 == 1:
                    left = 0.0 if i == 0 else u[i-1, j]
                    right = 0.0 if i == Nx-1 else u[i+1, j]
                    bottom = 0.0 if j == 0 else u[i, j-1]
                    top = 0.0 if j == Ny-1 else u[i, j+1]
                    
                    laplacian_u = (left + right - 2*u[i, j]) / hx2 + (bottom + top - 2*u[i, j]) / hy2
                    
                    u[i, j] += tau * (f_matrix[i, j] + laplacian_u)
        
        residual = np.max(np.abs(u - u_old))
        history.append(residual)
        
        if residual < tol:
            execution_time = time.time() - start_time
            print(f"✓ Метод сошелся за {iteration+1} итераций")
            print(f"Время: {execution_time:.4f} сек, невязка: {residual:.2e}")
            return u#, history, execution_time
        
        if iteration % 100 == 0:
            print(f"Итерация {iteration:4d}, невязка: {residual:.2e}")
            
        if np.isnan(residual) or residual > 1e10:
            print(f"✗ Метод расходится на итерации {iteration}")
            break
    
    execution_time = time.time() - start_time
    print(f"✗ Метод не сошелся за {max_iter} итераций")
    return u#, history, execution_time
# =============================================================================
# Шахматный метод (Red-Black ordering) для уравнения Пуассона
# =============================================================================
def checkerboard_method_poisson_solver(f_matrix, Nx, Ny, Lx=1.0, Ly=1.0, max_iter=1000, tol=1e-6):
    """
    Шахматный метод (Red-Black ordering) для уравнения Пуассона
    """
    start_time = time.time()

    # Раздельные шаги сетки
    total_points_x = Nx + 2
    total_points_y = Ny + 2
    hx = Lx / (total_points_x - 1)
    hy = Ly / (total_points_y - 1)
    hx2 = hx**2
    hy2 = hy**2
    
    print(f"Улучшенный шахматный метод:")
    print(f"Область: {Lx} x {Ly}")
    print(f"Внутренние точки: {Nx} x {Ny}")
    print(f"Шаги: hx = {hx:.4f}, hy = {hy:.4f}")
    print(f"hx² = {hx2:.6f}, hy² = {hy2:.6f}")
    
    u = np.zeros((Nx, Ny))
    history = []
    
    for iteration in range(max_iter):
        u_old = u.copy()
        
        # Черные узлы
        for i in range(Nx):
            for j in range(Ny):
                if (i + j) % 2 == 0:
                    left = 0.0 if i == 0 else u[i-1, j]
                    right = 0.0 if i == Nx-1 else u[i+1, j]
                    bottom = 0.0 if j == 0 else u[i, j-1]
                    top = 0.0 if j == Ny-1 else u[i, j+1]
                    
                    # Для прямоугольной области используем взвешенное среднее
                    u[i, j] = (hx2 * hy2 * f_matrix[i, j] + 
                               hy2 * (left + right) + 
                               hx2 * (bottom + top)) / (2 * (hx2 + hy2))
        
        # Красные узлы
        for i in range(Nx):
            for j in range(Ny):
                if (i + j) % 2 == 1:
                    left = 0.0 if i == 0 else u[i-1, j]
                    right = 0.0 if i == Nx-1 else u[i+1, j]
                    bottom = 0.0 if j == 0 else u[i, j-1]
                    top = 0.0 if j == Ny-1 else u[i, j+1]
                    
                    u[i, j] = (hx2 * hy2 * f_matrix[i, j] + 
                               hy2 * (left + right) + 
                               hx2 * (bottom + top)) / (2 * (hx2 + hy2))
        
        residual = np.max(np.abs(u - u_old))
        history.append(residual)
        
        if residual < tol:
            execution_time = time.time() - start_time
            print(f"✓ Метод сошелся за {iteration+1} итераций")
            print(f"Время: {execution_time:.4f} сек, невязка: {residual:.2e}")
            return u#, history, execution_time
        
        if iteration % 100 == 0:
            print(f"Итерация {iteration:4d}, невязка: {residual:.2e}")
            
        if np.isnan(residual) or residual > 1e10:
            print(f"✗ Метод расходится на итерации {iteration}")
            break
    
    execution_time = time.time() - start_time
    print(f"✗ Метод не сошелся за {max_iter} итераций")
    return u#, history, execution_time
# =============================================================================
# 
# =============================================================================

class CorrectPoissonSolver2D:
    """
    Решатель уравнения Пуассона с ТОЧНЫМИ формулами из файла
    """
    
    def __init__(self, Nx, Ny, Lx, Ly):
        self.Nx = Nx  # Количество интервалов по x
        self.Ny = Ny  # Количество интервалов по y
        self.Lx = Lx  # Длина  по x
        self.Ly = Ly  # Длина  по y
        self.hx = Lx / Nx #Шаг  по x
        self.hy = Ly / Ny  #Шаг  по y
        
        # Общее число точек включая границы
        self.Nx_total = Nx + 1  # 0, 1, 2, ..., Nx
        self.Ny_total = Ny + 1  # 0, 1, 2, ..., Ny        
        """
        Решение уравнения Пуассона Δu = -f
        """
    def adding_boundary_conditions(self,u):        
        u_f = np.zeros((self.N1+1, self.N2+1))
        for i in range(1,self.N1_total,1):
            for j in range(1,self.N2_total,1):
                u_f[i,j] = u[i-1,j-1]
        return u_f
    
    def solve_DST(self, f_inner):
        self.u_dst = solve_u_DST(f_inner, self.Nx, self.Ny,self.hx, self.hy, self.Lx, self.Ly)
        return self.u_dst

    def solve_ADI(self, f_inner):
        self.u_adi = adi_poisson_solver(f_inner, self.Nx-1, self.Ny-1, self.Lx, self.Ly,  1000, 1e-6)
        return self.u_adi
    
    def solve_triangular_method(self,f_inner):
        self.u_triangular_method = alternating_triangle_method_poisson_solver(f_inner, self.Nx-1, self.Ny-1, self.Lx, self.Ly,  1000, 1e-6)
        return self.u_triangular_method

    def solve_checkerboard_method(self,f_inner):
        self.u_checkerboard_method = checkerboard_method_poisson_solver(f_inner, self.Nx-1, self.Ny-1, self.Lx, self.Ly,  1000, 1e-6)
        return self.u_checkerboard_method
    
# Тестовая правая часть
def f_test(x, y):
    return -2*(y**2-y+x**2-x)#2 * np.pi**2 * np.sin(np.pi * x) * np.sin(np.pi * y)

def demonstrate_correct_formulas():
    """Демонстрация работы с ТОЧНЫМИ формулами из файла"""
    
    print("Решение уравнения Пуассона с ТОЧНЫМИ формулами из файла")
    print("=" * 60)
    
    # Параметры сетки
    N1, N2 = 16, 16  # Уменьшим для скорости демонстрации
    L1, L2 = 1.0, 1.0
        
    # Создаем решатель
    solver = CorrectPoissonSolver2D(N1, N2, L1, L2)
    
    print(f"Параметры сетки:")
    print(f"  N1 = {N1}, N2 = {N2}")
    print(f"  h1 = {solver.hx:.4f}, h2 = {solver.hy:.4f}")
    print(f"  Общее число точек: {solver.Nx_total} x {solver.Ny_total}")

    h1 = L1/N1
    h2 = L2/N2
    # Шаг 1: Дискретизация(Приближение/Сеточное задавание правой части)
    f = np.zeros((N1-1, N2-1))
    for i in range(0, N1-1):  # i = 1, 2, ..., N1-1
        for j in range(0, N2-1):  # j = 1, 2, ..., N2-1
            x = (i+1) * h1
            y = (j+1) * h2
            f[i, j] = f_test(x, y)
            #print(f[i, j], i,j)
    # Берем только внутренние точки
    #f_inner = f[1:N1, 1:N2]
    #print(len(f))  
    # Решаем уравнение
    print("\nРешаем уравнение Пуассона DST, ADI")
    solution_DST = solver.solve_DST(f)
    solution_ADI = solver.solve_ADI(f)
    #solution_triangular = solver.solve_triangular_method(f)
    #solution_checkerboard = solver.solve_checkerboard_method(f)
    
    # Аналитическое решение
    analytical = np.zeros((solver.Nx-1, solver.Ny-1))
    for i in range(0,solver.Nx-1):
        for j in range(0,solver.Ny-1):
            x = (i+1) * solver.hx
            y = (j+1) * solver.hy
            analytical[i, j] = (x**2-x)*(y**2-y)#np.sin(np.pi * x) * np.sin(np.pi * y)
    #print(analytical)
            
    # Вычисляем погрешность
    error_DST = np.abs(solution_DST - analytical)    
    print(f"\nРезультаты DST:")
    print(f"Максимальная погрешность: {np.max(error_DST):.2e}")
    print(f"Средняя погрешность: {np.mean(error_DST):.2e}")    
    # Визуализация
    visualize_solutions(solution_DST, analytical, error_DST, L1, L2)
    
    # Вычисляем погрешность
    error_ADI = np.abs(solution_ADI - analytical)    
    print(f"\nРезультаты ADI:")
    print(f"Максимальная погрешность: {np.max(error_ADI):.2e}")
    print(f"Средняя погрешность: {np.mean(error_ADI):.2e}")    
    # Визуализация
    visualize_solutions(solution_ADI, analytical, error_ADI, L1, L2)
    
    error = np.abs(solution_ADI - solution_DST)
    print(f"\nРезультаты Разницы решений:")
    print(f"Максимальная погрешность: {np.max(error):.2e}")
    print(f"Средняя погрешность: {np.mean(error):.2e}")  
    '''

    # Вычисляем погрешность
    error_triangular = np.abs(solution_triangular - analytical)    
    print(f"\nРезультаты Попеременно-треугольный метод:")
    print(f"Максимальная погрешность: {np.max(error_triangular):.2e}")
    print(f"Средняя погрешность: {np.mean(error_triangular):.2e}")    
    # Визуализация
    visualize_solutions(solution_triangular, analytical, error_triangular, L1, L2)
    
    error = np.abs(solution_triangular - solution_DST)
    print(f"\nРезультаты Разницы решений DST - Попеременно-треугольный :")
    print(f"Максимальная погрешность: {np.max(error):.2e}")
    print(f"Средняя погрешность: {np.mean(error):.2e}")
    
    # Вычисляем погрешность
    error_checkerboard = np.abs(solution_checkerboard - analytical)    
    print(f"\nРезультаты Шахматный метод:")
    print(f"Максимальная погрешность: {np.max(error_checkerboard):.2e}")
    print(f"Средняя погрешность: {np.mean(error_checkerboard):.2e}")    
    # Визуализация
    visualize_solutions(solution_checkerboard, analytical, error_checkerboard, L1, L2)

    error = np.abs(solution_checkerboard - solution_DST)
    print(f"\nРезультаты Разницы решений DST - Шахматный метод:")
    print(f"Максимальная погрешность: {np.max(error):.2e}")
    print(f"Средняя погрешность: {np.mean(error):.2e}")  
    '''
    
# Запуск демонстрации
if __name__ == "__main__":
    
    demonstrate_correct_formulas()    
    #check_lamda()









