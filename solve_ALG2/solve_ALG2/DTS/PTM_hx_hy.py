import numpy as np
import matplotlib.pyplot as plt
import time

class PoissonSolverPTM:
    def __init__(self, Nx=11, Ny=11, Lx=1.0, Ly=1.0, tau_factor=0.5):
        """
        Nx, Ny - количество узлов (включая границы)
        Lx, Ly - размеры области
        tau_factor - коэффициент для tau (итерационного параметра)
        """
        self.Nx = Nx  # всего точек по x, включая границы
        self.Ny = Ny  # всего точек по y, включая границы
        self.Lx = Lx
        self.Ly = Ly
        
        # Сетка
        self.hx = Lx / (Nx) if Nx > 1 else Lx
        self.hy = Ly / (Ny) if Ny > 1 else Ly
        
        # Параметры метода
        self.tau = 0.01#tau_factor * self.hx * self.hy  # итерационный параметр
        
        self.hx2 = self.hx * self.hx
        self.hy2 = self.hy * self.hy
        self.r = self.hx2 / self.hy2 if self.hy2 != 0 else 1.0
        self.q = 1.0 / self.r if self.r != 0 else 1.0
        self.ax = 2*self.hx2 / self.tau if self.tau != 0 else 1.0
        self.ay = 2*self.hy2 / self.tau if self.tau != 0 else 1.0
        self.o_ax2 = 1.0 / (self.ax + 2)
        self.o_ay2 = 1.0 / (self.ay + 2)
        
        # Массивы
        self.u = np.zeros((Nx, Ny))
        self.f = np.zeros((Nx, Ny))
        self.u_exact = np.zeros((Nx, Ny))
        
        # Для ускорения вычислений
        self.Fmax = 0.0
        
    def set_test_problem(self, test_case=1):
        """Тестовые задачи"""
        x = np.linspace(0, self.Lx, self.Nx)
        y = np.linspace(0, self.Ly, self.Ny)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        if test_case == 1:
            # Тест 1: u = sin(πx)sin(πy)
            # Уравнение: Δu = -2π² sin(πx)sin(πy) = f
            self.u_exact = np.sin(np.pi * X) * np.sin(np.pi * Y)
            self.f = -2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)  # Δu = f
            self.problem_name = "u = sin(πx)sin(πy)"
            
        elif test_case == 2:
            # Тест 2: u = x(1-x)y(1-y)
            # Уравнение: Δu = -2[x(1-x) + y(1-y)] = f
            self.u_exact = X * (1 - X) * Y * (1 - Y)
            self.f = -2.0 * (X * (1 - X) + Y * (1 - Y))  # Δu = f
            self.problem_name = "u = x(1-x)y(1-y)"
            
        elif test_case == 3:
            # Тест 3: u = (x² - x)(y² - y)
            # Δu = 2(y² - y) + 2(x² - x) = 2(x² - x + y² - y)
            # f = Δu = 2(x² - x + y² - y)
            self.u_exact = (X*X - X) * (Y*Y - Y)
            self.f = 2.0 * (X*X - X + Y*Y - Y)  # Δu = f
            self.problem_name = "u = (x²-x)(y²-y)"
            
        elif test_case == 4:
            # Тест 4: простая константа для проверки
            self.u_exact = np.ones_like(X)
            self.f = np.zeros_like(X)
            self.problem_name = "u = 1 (константа)"
            
        # Начальное приближение (нули)
        self.u = np.zeros_like(self.u_exact)
        
        # Граничные условия (Дирихле) - нули для тестов 1-3
        if test_case in [1, 2, 3]:
            self.u[0, :] = 0.0
            self.u[-1, :] = 0.0
            self.u[:, 0] = 0.0
            self.u[:, -1] = 0.0
        else:
            self.u[0, :] = self.u_exact[0, :]
            self.u[-1, :] = self.u_exact[-1, :]
            self.u[:, 0] = self.u_exact[:, 0]
            self.u[:, -1] = self.u_exact[:, -1]
        
        # Вычисляем Fmax = max|f * hx²|
        self.Fmax = np.max(np.abs(self.f * self.hx2))
        if self.Fmax < 1e-15:
            self.Fmax = 1.0  # чтобы избежать деления на ноль
            
    def compute_max_residual(self):
        """Вычисление максимальной невязки: max|L(u) - f|"""
        max_res = 0.0
        
        for i in range(1, self.Nx-1):
            for j in range(1, self.Ny-1):
                # L(u) = (u_{i+1,j} + u_{i-1,j} - 2u_{ij})/hx² + (u_{i,j+1} + u_{i,j-1} - 2u_{ij})/hy²
                # Умножаем на hx²: (u_{i+1,j} + u_{i-1,j} - 2u_{ij}) + r*(u_{i,j+1} + u_{i,j-1} - 2u_{ij})
                residual = (self.u[i+1, j] + self.u[i-1, j] - 2*self.u[i, j]) + \
                          self.r * (self.u[i, j+1] + self.u[i, j-1] - 2*self.u[i, j]) - \
                          self.f[i, j] * self.hx2
                
                abs_res = abs(residual)
                if abs_res > max_res:
                    max_res = abs_res
        
        return max_res
    
    def compute_error(self):
        """Вычисление максимальной ошибки относительно точного решения"""
        max_err = 0.0
        for i in range(self.Nx):
            for j in range(self.Ny):
                err = abs(self.u[i, j] - self.u_exact[i, j])
                if err > max_err:
                    max_err = err
        return max_err
    
    def solve(self, max_iter=100000, EPS=1e-6, verbose=True):
        """Решение ПТМ методом с критерием остановки: Error0 < EPS * Fmax"""
        start_time = time.time()
        
        # Сохраняем граничные значения
        u_boundaries = np.zeros((self.Nx, self.Ny))
        u_boundaries[0, :] = self.u[0, :].copy()
        u_boundaries[-1, :] = self.u[-1, :].copy()
        u_boundaries[:, 0] = self.u[:, 0].copy()
        u_boundaries[:, -1] = self.u[:, -1].copy()
        
        print(f"Начало решения: Nx={self.Nx}, Ny={self.Ny}, hx={self.hx:.4f}, hy={self.hy:.4f}")
        print(f"Критерий остановки: max_residual < {EPS} * {self.Fmax:.6e} = {EPS * self.Fmax:.6e}")
        
        errors = []
        residuals = []
        
        for k in range(max_iter):
            u_old = self.u.copy()
            
            # Шаг 1: слева направо, снизу вверх
            for j in range(1, self.Ny-1):
                for i in range(1, self.Nx-1):
                    self.u[i, j] = (self.u[i-1, j] + self.u[i, j-1] + 
                                   self.ax * u_old[i, j] + 
                                   self.r * (u_old[i, j+1] + u_old[i, j-1] - 2*u_old[i, j]) - 
                                   self.f[i, j] * self.hx2) * self.o_ax2
            
            u_old = self.u.copy()  # Теперь u_old = u^{n+1/2}
            
            # Шаг 2: справа налево, сверху вниз
            for j in range(self.Ny-2, 0, -1):
                for i in range(self.Nx-2, 0, -1):
                    self.u[i, j] = (self.u[i+1, j] + self.u[i, j+1] + 
                                   self.ay * u_old[i, j] + 
                                   self.q * (u_old[i+1, j] + u_old[i-1, j] - 2*u_old[i, j]) - 
                                   self.f[i, j] * self.hy2) * self.o_ay2
            
            # Восстанавливаем граничные условия
            self.u[0, :] = u_boundaries[0, :]
            self.u[-1, :] = u_boundaries[-1, :]
            self.u[:, 0] = u_boundaries[:, 0]
            self.u[:, -1] = u_boundaries[:, -1]
            
            # Вычисление максимальной невязки (Error0)
            Error0 = self.compute_max_residual()
            
            # Вычисление ошибки
            error = self.compute_error()
            
            errors.append(error)
            residuals.append(Error0)
            
            if verbose and k % 1000 == 0:
                print(f"Итерация {k:6d}: max_residual = {Error0:.3e}, error = {error:.3e}")
            
            # Критерий остановки: Error0 < EPS * Fmax
            if Error0 < EPS * self.Fmax:
                elapsed_time = time.time() - start_time
                print(f"\n✓ Сходимость достигнута за {k+1} итераций ({elapsed_time:.2f} секунд)")
                print(f"  Финальная невязка = {Error0:.3e}")
                print(f"  Финальная ошибка = {error:.3e}")
                break
                
            if k == max_iter - 1:
                elapsed_time = time.time() - start_time
                print(f"\n✗ Достигнуто максимальное число итераций ({max_iter}) за {elapsed_time:.2f} секунд")
                print(f"  Финальная невязка = {Error0:.3e}")
                print(f"  Финальная ошибка = {error:.3e}")
        
        self.errors = errors
        self.residuals = residuals
        self.iterations = k + 1
        self.EPS = EPS
        
        return errors, residuals
    
    def plot_solution(self):
        """Графики решения"""
        fig = plt.figure(figsize=(15, 10))
        
        # Создаем координатные сетки
        x = np.linspace(0, self.Lx, self.Nx)
        y = np.linspace(0, self.Ly, self.Ny)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        # 1. Точное решение
        ax1 = plt.subplot(2, 3, 1)
        contour1 = ax1.contourf(X, Y, self.u_exact, 20, cmap='viridis')
        ax1.set_title(f'Точное решение\n{self.problem_name}')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_aspect('equal')
        plt.colorbar(contour1, ax=ax1)
        
        # 2. Численное решение
        ax2 = plt.subplot(2, 3, 2)
        contour2 = ax2.contourf(X, Y, self.u, 20, cmap='viridis')
        ax2.set_title(f'Численное решение\nИтераций: {self.iterations}')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_aspect('equal')
        plt.colorbar(contour2, ax=ax2)
        
        # 3. Ошибка (разность)
        ax3 = plt.subplot(2, 3, 3)
        error = self.u - self.u_exact
        max_abs_error = np.max(np.abs(error))
        contour3 = ax3.contourf(X, Y, error, 20, cmap='RdBu', vmin=-max_abs_error, vmax=max_abs_error)
        ax3.set_title(f'Ошибка (числ. - точн.)\nМакс. ошибка: {max_abs_error:.2e}')
        ax3.set_xlabel('x')
        ax3.set_ylabel('y')
        ax3.set_aspect('equal')
        plt.colorbar(contour3, ax=ax3)
        
        # 4. Правая часть f(x,y)
        ax4 = plt.subplot(2, 3, 4)
        contour4 = ax4.contourf(X, Y, self.f, 20, cmap='plasma')
        ax4.set_title('Правая часть f(x,y)')
        ax4.set_xlabel('x')
        ax4.set_ylabel('y')
        ax4.set_aspect('equal')
        plt.colorbar(contour4, ax=ax4)
        
        # 5. Сходимость невязки
        ax5 = plt.subplot(2, 3, 5)
        iterations = range(len(self.residuals))
        ax5.semilogy(iterations, self.residuals, 'b-', linewidth=2, label='Невязка')
        ax5.axhline(y=self.EPS * self.Fmax, color='r', linestyle='--', 
                   label=f'Критерий: {self.EPS}*Fmax = {self.EPS*self.Fmax:.2e}')
        ax5.set_title('Сходимость невязки')
        ax5.set_xlabel('Итерация')
        ax5.set_ylabel('Макс. невязка (log scale)')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
        
        # 6. Сходимость ошибки
        ax6 = plt.subplot(2, 3, 6)
        ax6.semilogy(iterations, self.errors, 'g-', linewidth=2)
        ax6.set_title('Сходимость ошибки')
        ax6.set_xlabel('Итерация')
        ax6.set_ylabel('Макс. ошибка (log scale)')
        ax6.grid(True, alpha=0.3)
        
        plt.suptitle(f'ПТМ метод: {self.problem_name}, Сетка: {self.Nx}×{self.Ny}, hx={self.hx:.4f}, hy={self.hy:.4f}', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.show()
    
    def print_statistics(self):
        """Вывод статистики"""
        print("\n" + "="*60)
        print("РЕЗУЛЬТАТЫ РАСЧЕТА:")
        print("="*60)
        print(f"Задача: {self.problem_name}")
        print(f"Область: {self.Lx}×{self.Ly}")
        print(f"Сетка: {self.Nx}×{self.Ny} точек")
        print(f"Шаги сетки: hx = {self.hx:.4f}, hy = {self.hy:.4f}")
        print(f"Параметр τ: {self.tau:.6f}")
        print(f"Коэффициенты: ax={self.ax:.2f}, ay={self.ay:.2f}, r={self.r:.2f}")
        print(f"Выполнено итераций: {self.iterations}")
        
        max_err = np.max(np.abs(self.u - self.u_exact))
        mean_err = np.mean(np.abs(self.u - self.u_exact))
        residual_norm = self.compute_max_residual()
        
        print(f"\nПогрешности:")
        print(f"  Максимальная погрешность: {max_err:.6e}")
        print(f"  Средняя погрешность: {mean_err:.6e}")
        print(f"  Норма невязки: {residual_norm:.6e}")
        
        # Проверка критерия
        criterion = self.EPS * self.Fmax
        print(f"\nКритерий остановки (EPS={self.EPS}):")
        print(f"  Fmax = max|f*hx²| = {self.Fmax:.6e}")
        print(f"  EPS * Fmax = {criterion:.6e}")
        print(f"  Фактическая невязка = {residual_norm:.6e}")
        print(f"  Условие выполнено: {'ДА' if residual_norm < criterion else 'НЕТ'}")
        print("="*60)

# Константы
EPS = 1e-6  # Требуемая точность

# Тестируем
if __name__ == "__main__":
    print("=" * 60)
    print("ТЕСТ 1: Квадратная область 1×1")
    print("Задача: u = sin(πx)sin(πy)")
    print("Сетка: 11×11 точек")
    print("=" * 60)
    
    # Квадратная область: 11×11 точек (быстро)
    solver1 = PoissonSolverPTM(Nx=11, Ny=11, Lx=1.0, Ly=1.0, tau_factor=0.5)
    solver1.set_test_problem(test_case=1)
    
    errors1, residuals1 = solver1.solve(max_iter=10000, EPS=EPS, verbose=True)
    solver1.print_statistics()
    solver1.plot_solution()
    
    print("\n\n" + "=" * 60)
    print("ТЕСТ 2: Прямоугольная область 2×1")
    print("Задача: u = sin(πx)sin(πy)")
    print("Сетка: 21×11 точек")
    print("=" * 60)
    
    # Прямоугольная область: 21×11 точек
    solver2 = PoissonSolverPTM(Nx=21, Ny=11, Lx=2.0, Ly=1.0, tau_factor=0.5)
    solver2.set_test_problem(test_case=1)
    
    errors2, residuals2 = solver2.solve(max_iter=10000, EPS=EPS, verbose=True)
    solver2.print_statistics()
    solver2.plot_solution()
    
    print("\n\n" + "=" * 60)
    print("ТЕСТ 3: Квадратная область 1×1")
    print("Задача: u = (x² - x)(y² - y)")
    print("Сетка: 11×11 точек")
    print("=" * 60)
    
    solver3 = PoissonSolverPTM(Nx=11, Ny=11, Lx=1.0, Ly=1.0, tau_factor=0.5)
    solver3.set_test_problem(test_case=3)
    
    errors3, residuals3 = solver3.solve(max_iter=10000, EPS=EPS, verbose=True)
    solver3.print_statistics()
    solver3.plot_solution()
