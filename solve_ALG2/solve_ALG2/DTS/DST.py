import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def lamda_k_0(k,h,N):
    
    return 4*(np.sin(k*np.pi/(2*N))**2)/(h**2)

def phi_k2(f_ij,i,k2,N2):
    SUM = 0.0
   # print(len(f_ij[i]))
    for j in range(1,N2,1):
        #print(j)
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
    """Правильная формула для собственных значений"""
    return 4/(h**2) * np.sin(np.pi * k * h / (2 * L))**2

def solve_u(f_inner, N1, N2, h1, h2, L1, L2):
    
    u = np.zeros((N1-1, N2-1))
    print(N1-1,N2-1)
    for i in range(1, N1):
        for j in range(1, N2):
            u[i-1,j-1] = u_ij(f_inner,i,j, h1,h2, N1,N2)*4/(N1*N2)           
    return u
      
class CorrectPoissonSolver2D:
    """
    Решатель уравнения Пуассона с ТОЧНЫМИ формулами из файла
    """
    
    def __init__(self, N1, N2, L1, L2):
        self.N1 = N1  # Количество интервалов по x
        self.N2 = N2  # Количество интервалов по y
        self.L1 = L1  # Длина  по x
        self.L2 = L2  # Длина  по y
        self.h1 = L1 / N1 #Шаг  по x
        self.h2 = L2 / N2  #Шаг  по x
        
        # Общее число точек включая границы
        self.N1_total = N1 + 1  # 0, 1, 2, ..., N1
        self.N2_total = N2 + 1  # 0, 1, 2, ..., N2
    
    def solve(self, f_inner):
        """
        Решение уравнения Пуассона Δu = -f с использованием ТОЧНЫХ формул
        """
        #print('f=',f_inner)
        u = solve_u(f_inner, self.N1, self.N2,self.h1, self.h2, self.L1, self.L2)
        #print('u =', u)           
        '''
        u_f = np.zeros((self.N1+1, self.N2+1))
        for i in range(1,self.N1_total,1):
            for j in range(1,self.N2_total,1):
                u_f[i,j] = u[i-1,j-1]
        #print(u_)
        '''
        return u
    
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
    print(f"  h1 = {solver.h1:.4f}, h2 = {solver.h2:.4f}")
    print(f"  Общее число точек: {solver.N1_total} x {solver.N2_total}")

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
    print(len(f))  
    # Решаем уравнение
    print("\nРешаем уравнение Пуассона...")
    solution = solver.solve(f)
    
    # Аналитическое решение
    analytical = np.zeros((solver.N1-1, solver.N2-1))
    for i in range(0,solver.N1-1):
        for j in range(0,solver.N2-1):
            x = (i+1) * solver.h1
            y = (j+1) * solver.h2
            analytical[i, j] = (x**2-x)*(y**2-y)#np.sin(np.pi * x) * np.sin(np.pi * y)
    print(analytical)
    # Вычисляем погрешность
    error = np.abs(solution - analytical)
    
    print(f"\nРезультаты:")
    print(f"Максимальная погрешность: {np.max(error):.2e}")
    print(f"Средняя погрешность: {np.mean(error):.2e}")
    
    # Визуализация
    visualize_solutions(solution, analytical, error, L1, L2)
    
    return solution, analytical, error

def visualize_solutions(numerical, analytical, error, L1, L2):
    """Визуализация численного и аналитического решений"""
    
    N1_total, N2_total = numerical.shape
    x = np.linspace(0, L1, N1_total)
    y = np.linspace(0, L2, N2_total)
    X, Y = np.meshgrid(x, y, indexing='ij')
    
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
    '''
    # Аналитическое решение
    ax3 = fig.add_subplot(133, projection='3d')
    surf3 = ax3.plot_surface(X, Y, error, cmap='hot', alpha=0.8)
    ax3.set_title('F')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('f')
    '''
    # Погрешность
    ax3 = fig.add_subplot(133, projection='3d')
    surf3 = ax3.plot_surface(X, Y, error, cmap='hot', alpha=0.8)
    ax3.set_title('Погрешность')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    ax3.set_zlabel('|error|')
    
    plt.tight_layout()
    plt.show()
    
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
    #contour() - только линии контуров
    #contourf() - залитые области между контурами

    
# Запуск демонстрации
if __name__ == "__main__":
    
    demonstrate_correct_formulas()    
    #check_lamda()









