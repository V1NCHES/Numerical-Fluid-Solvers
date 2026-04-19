import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def load_data_from_file(filename):
    """Загружает данные из файла, где каждая строка содержит x, y, значение"""
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line:  # пропускаем пустые строки
                parts = line.split()
                if len(parts) == 3:
                    x, y, z = map(float, parts)
                    data.append((x, y, z))
    
    # Преобразуем в numpy массивы
    x_vals = np.array([d[0] for d in data])
    y_vals = np.array([d[1] for d in data])
    z_vals = np.array([d[2] for d in data])
    
    # Находим уникальные значения x и y
    unique_x = np.unique(x_vals)
    unique_y = np.unique(y_vals)
    
    # Создаем сетку для тепловой карты
    X, Y = np.meshgrid(unique_x, unique_y)
    Z = np.zeros_like(X)
    
    # Заполняем матрицу Z значениями
    for x, y, z in data:
        i = np.where(unique_x == x)[0][0]
        j = np.where(unique_y == y)[0][0]
        Z[j, i] = z  # j - индекс строки, i - индекс столбца
    
    return X, Y, Z

def plot_comparison_with_difference(file1, file2, title1="Аналитическое решение", title2="Треугольное решение"):
    """Строит тепловые карты, контуры и разность для двух файлов"""
    
    # Загружаем данные
    X1, Y1, Z1 = load_data_from_file(file1)
    X2, Y2, Z2 = load_data_from_file(file2)
    
    # Вычисляем разность (Z1 - Z2)
    Z_diff = Z1 - Z2
    
    # Вычисляем относительную разность в процентах
    # Добавляем небольшое значение для избежания деления на ноль
    Z_rel_diff = np.abs(Z_diff) / (np.abs(Z1) + 1e-10) * 100
    
    # Создаем фигуру с 6 подграфиками (2x3)
    fig = plt.figure(figsize=(20, 14))
    
    # 1. Тепловая карта для файла 1
    ax1 = fig.add_subplot(2, 3, 1)
    contour1 = ax1.contourf(X1, Y1, Z1, levels=50, cmap='viridis')
    ax1.set_title(f"{title1}\nТепловая карта", fontsize=12)
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    plt.colorbar(contour1, ax=ax1, label='Значение')
    ax1.set_aspect('equal')
    
    # 2. Тепловая карта для файла 2
    ax2 = fig.add_subplot(2, 3, 2)
    contour2 = ax2.contourf(X2, Y2, Z2, levels=50, cmap='viridis')
    ax2.set_title(f"{title2}\nТепловая карта", fontsize=12)
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    plt.colorbar(contour2, ax=ax2, label='Значение')
    ax2.set_aspect('equal')
    
    # 3. Абсолютная разность
    ax3 = fig.add_subplot(2, 3, 3)
    # Используем diverging colormap для разности
    max_abs_diff = np.max(np.abs(Z_diff))
    levels = np.linspace(-max_abs_diff, max_abs_diff, 50)
    contour3 = ax3.contourf(X1, Y1, Z_diff, levels=levels, cmap='RdBu_r')
    ax3.set_title("Абсолютная разность\n(Z1 - Z2)", fontsize=12)
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    plt.colorbar(contour3, ax=ax3, label='Разность')
    ax3.set_aspect('equal')
    
    # 4. Контуры для файла 1
    ax4 = fig.add_subplot(2, 3, 4)
    contour4 = ax4.contour(X1, Y1, Z1, levels=15, colors='black', linewidths=0.8)
    ax4.clabel(contour4, inline=True, fontsize=8)
    ax4.set_title(f"{title1}\nКонтуры", fontsize=12)
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_aspect('equal')
    
    # 5. Контуры для файла 2
    ax5 = fig.add_subplot(2, 3, 5)
    contour5 = ax5.contour(X2, Y2, Z2, levels=15, colors='black', linewidths=0.8)
    ax5.clabel(contour5, inline=True, fontsize=8)
    ax5.set_title(f"{title2}\nКонтуры", fontsize=12)
    ax5.set_xlabel('X')
    ax5.set_ylabel('Y')
    ax5.set_aspect('equal')
    
    # 6. Относительная разность в процентах
    ax6 = fig.add_subplot(2, 3, 6)
    contour6 = ax6.contourf(X1, Y1, Z_rel_diff, levels=50, cmap='hot')
    ax6.set_title("Относительная разность\n|Z1 - Z2| / |Z1| × 100%", fontsize=12)
    ax6.set_xlabel('X')
    ax6.set_ylabel('Y')
    plt.colorbar(contour6, ax=ax6, label='Разность (%)')
    ax6.set_aspect('equal')
    
    plt.tight_layout()
    plt.show()
    
    # Отдельная фигура с подробной статистикой разности
    fig2, ((ax7, ax8), (ax9, ax10)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # 7. Абсолютная разность с контурами
    im7 = ax7.contourf(X1, Y1, Z_diff, levels=50, cmap='RdBu_r')
    contour7 = ax7.contour(X1, Y1, Z_diff, levels=10, colors='black', linewidths=0.5, alpha=0.7)
    ax7.clabel(contour7, inline=True, fontsize=8)
    ax7.set_title("Абсолютная разность с контурами", fontsize=12)
    ax7.set_xlabel('X')
    ax7.set_ylabel('Y')
    ax7.set_aspect('equal')
    plt.colorbar(im7, ax=ax7, label='Z1 - Z2')
    
    # 8. Относительная разность с контурами
    im8 = ax8.contourf(X1, Y1, Z_rel_diff, levels=50, cmap='YlOrRd')
    contour8 = ax8.contour(X1, Y1, Z_rel_diff, levels=10, colors='black', linewidths=0.5, alpha=0.7)
    ax8.clabel(contour8, inline=True, fontsize=8)
    ax8.set_title("Относительная разность (%) с контурами", fontsize=12)
    ax8.set_xlabel('X')
    ax8.set_ylabel('Y')
    ax8.set_aspect('equal')
    plt.colorbar(im8, ax=ax8, label='Разность (%)')
    
    # 9. Гистограмма абсолютной разности
    ax9.hist(Z_diff.flatten(), bins=50, edgecolor='black', alpha=0.7)
    ax9.axvline(x=0, color='r', linestyle='--', linewidth=1, label='Нулевая разность')
    ax9.set_title("Распределение абсолютной разности", fontsize=12)
    ax9.set_xlabel('Z1 - Z2')
    ax9.set_ylabel('Частота')
    ax9.legend()
    ax9.grid(True, alpha=0.3)
    
    # 10. Гистограмма относительной разности
    ax10.hist(Z_rel_diff.flatten(), bins=50, edgecolor='black', alpha=0.7, color='orange')
    ax10.set_title("Распределение относительной разности (%)", fontsize=12)
    ax10.set_xlabel('Относительная разность (%)')
    ax10.set_ylabel('Частота')
    ax10.grid(True, alpha=0.3)
    
    # Вывод статистики
    print("=" * 60)
    print("СТАТИСТИКА СРАВНЕНИЯ")
    print("=" * 60)
    print(f"Максимальная абсолютная разность: {np.max(np.abs(Z_diff)):.6e}")
    print(f"Средняя абсолютная разность: {np.mean(np.abs(Z_diff)):.6e}")
    print(f"Средняя относительная разность: {np.mean(Z_rel_diff):.2f}%")
    print(f"Максимальная относительная разность: {np.max(Z_rel_diff):.2f}%")
    print(f"RMSE (среднеквадратичная ошибка): {np.sqrt(np.mean(Z_diff**2)):.6e}")
    print(f"Коэффициент корреляции: {np.corrcoef(Z1.flatten(), Z2.flatten())[0,1]:.6f}")
    print("=" * 60)
    
    plt.tight_layout()
    plt.show()
    
    # Третья фигура: 3D поверхности для визуализации разности
    from mpl_toolkits.mplot3d import Axes3D
    
    fig3 = plt.figure(figsize=(16, 6))
    
    # 3D поверхность Z1
    ax11 = fig3.add_subplot(1, 3, 1, projection='3d')
    surf1 = ax11.plot_surface(X1, Y1, Z1, cmap='viridis', alpha=0.8, edgecolor='none')
    ax11.set_title(f"{title1} - 3D", fontsize=12)
    ax11.set_xlabel('X')
    ax11.set_ylabel('Y')
    ax11.set_zlabel('Z')
    fig3.colorbar(surf1, ax=ax11, shrink=0.5, label='Значение')
    
    # 3D поверхность Z2
    ax12 = fig3.add_subplot(1, 3, 2, projection='3d')
    surf2 = ax12.plot_surface(X2, Y2, Z2, cmap='viridis', alpha=0.8, edgecolor='none')
    ax12.set_title(f"{title2} - 3D", fontsize=12)
    ax12.set_xlabel('X')
    ax12.set_ylabel('Y')
    ax12.set_zlabel('Z')
    fig3.colorbar(surf2, ax=ax12, shrink=0.5, label='Значение')
    
    # 3D поверхность разности
    ax13 = fig3.add_subplot(1, 3, 3, projection='3d')
    surf3 = ax13.plot_surface(X1, Y1, Z_diff, cmap='coolwarm', alpha=0.8, edgecolor='none')
    ax13.set_title("Разность (3D)", fontsize=12)
    ax13.set_xlabel('X')
    ax13.set_ylabel('Y')
    ax13.set_zlabel('Z1 - Z2')
    fig3.colorbar(surf3, ax=ax13, shrink=0.5, label='Разность')
    
    plt.tight_layout()
    plt.show()

# Использование функции
if __name__ == "__main__":
    # Укажите пути к вашим файлам
    file1 = "analytical_solution.txt"
    file2 = "solution_triangular.txt"
    
    # Построение графиков с разностью
    plot_comparison_with_difference(file1, file2)