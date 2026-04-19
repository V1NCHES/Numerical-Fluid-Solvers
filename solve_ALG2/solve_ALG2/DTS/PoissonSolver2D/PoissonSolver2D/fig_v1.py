import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

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
    
    return X, Y, Z, unique_x, unique_y

def plot_comparison_with_difference(file1, file2, title1="Аналитическое решение", title2="Треугольное решение"):
    """Строит тепловые карты, контуры и разность для двух файлов"""
    
    # Загружаем данные
    X1, Y1, Z1, ux1, uy1 = load_data_from_file(file1)
    X2, Y2, Z2, ux2, uy2 = load_data_from_file(file2)
    
    # Проверяем, что сетки совпадают
    if not np.allclose(ux1, ux2) or not np.allclose(uy1, uy2):
        print("Предупреждение: сетки координат в файлах не совпадают!")
    
    # Вычисляем разность (Z1 - Z2)
    Z_diff = Z1 - Z2
    
    # Вычисляем относительную разность в процентах
    # Добавляем небольшое значение для избежания деления на ноль
    Z_rel_diff = np.zeros_like(Z_diff)
    for i in range(Z_diff.shape[0]):
        for j in range(Z_diff.shape[1]):
            if abs(Z1[i, j]) > 1e-12:  # если значение не слишком близко к нулю
                Z_rel_diff[i, j] = abs(Z_diff[i, j]) / abs(Z1[i, j]) * 100
            elif abs(Z_diff[i, j]) > 1e-12:  # если разность значима, а Z1 ~ 0
                Z_rel_diff[i, j] = np.inf
            else:  # если и Z1 и разность близки к нулю
                Z_rel_diff[i, j] = 0.0
    
    # Создаем фигуру с 6 подграфиками (2x3)
    fig = plt.figure(figsize=(20, 14))
    
    # 1. Тепловая карта для файла 1
    ax1 = fig.add_subplot(2, 3, 1)
    contour1 = ax1.contourf(X1, Y1, Z1, levels=50, cmap='viridis')
    ax1.set_title(f"{title1}\nТепловая карта", fontsize=12, fontweight='bold')
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    plt.colorbar(contour1, ax=ax1, label='Значение')
    ax1.set_aspect('equal')
    # Добавляем точки для граничных условий
    ax1.scatter([ux1[0], ux1[-1]], [uy1[0], uy1[-1]], color='red', s=50, marker='o', 
                alpha=0.7, label='Граничные точки', zorder=5)
    ax1.legend(loc='upper right', fontsize=9)
    
    # 2. Тепловая карта для файла 2
    ax2 = fig.add_subplot(2, 3, 2)
    contour2 = ax2.contourf(X2, Y2, Z2, levels=50, cmap='viridis')
    ax2.set_title(f"{title2}\nТепловая карта", fontsize=12, fontweight='bold')
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    plt.colorbar(contour2, ax=ax2, label='Значение')
    ax2.set_aspect('equal')
    ax2.scatter([ux2[0], ux2[-1]], [uy2[0], uy2[-1]], color='red', s=50, marker='o', 
                alpha=0.7, label='Граничные точки', zorder=5)
    ax2.legend(loc='upper right', fontsize=9)
    
    # 3. Абсолютная разность
    ax3 = fig.add_subplot(2, 3, 3)
    # Используем diverging colormap для разности
    max_abs_diff = np.max(np.abs(Z_diff))
    if max_abs_diff > 0:
        levels = np.linspace(-max_abs_diff, max_abs_diff, 50)
        contour3 = ax3.contourf(X1, Y1, Z_diff, levels=levels, cmap='RdBu_r')
    else:
        contour3 = ax3.contourf(X1, Y1, Z_diff, cmap='RdBu_r')
    ax3.set_title("Абсолютная разность\n(Z1 - Z2)", fontsize=12, fontweight='bold')
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    plt.colorbar(contour3, ax=ax3, label='Разность')
    ax3.set_aspect('equal')
    
    # 4. Контуры для файла 1
    ax4 = fig.add_subplot(2, 3, 4)
    contour4 = ax4.contour(X1, Y1, Z1, levels=15, colors='black', linewidths=0.8)
    ax4.clabel(contour4, inline=True, fontsize=8)
    ax4.set_title(f"{title1}\nКонтуры", fontsize=12, fontweight='bold')
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_aspect('equal')
    # Заливка для контуров
    ax4.contourf(X1, Y1, Z1, levels=50, cmap='viridis', alpha=0.5)
    
    # 5. Контуры для файла 2
    ax5 = fig.add_subplot(2, 3, 5)
    contour5 = ax5.contour(X2, Y2, Z2, levels=15, colors='black', linewidths=0.8)
    ax5.clabel(contour5, inline=True, fontsize=8)
    ax5.set_title(f"{title2}\nКонтуры", fontsize=12, fontweight='bold')
    ax5.set_xlabel('X')
    ax5.set_ylabel('Y')
    ax5.set_aspect('equal')
    ax5.contourf(X2, Y2, Z2, levels=50, cmap='viridis', alpha=0.5)
    
    # 6. Относительная разность в процентах
    ax6 = fig.add_subplot(2, 3, 6)
    # Заменяем бесконечные значения на максимальные конечные
    Z_rel_diff_plot = Z_rel_diff.copy()
    finite_vals = Z_rel_diff_plot[np.isfinite(Z_rel_diff_plot)]
    if len(finite_vals) > 0:
        max_finite = np.max(finite_vals)
        Z_rel_diff_plot[~np.isfinite(Z_rel_diff_plot)] = max_finite * 1.1
    
    contour6 = ax6.contourf(X1, Y1, Z_rel_diff_plot, levels=50, cmap='hot')
    ax6.set_title("Относительная разность\n|Z1 - Z2| / |Z1| × 100%", fontsize=12, fontweight='bold')
    ax6.set_xlabel('X')
    ax6.set_ylabel('Y')
    plt.colorbar(contour6, ax=ax6, label='Разность (%)')
    ax6.set_aspect('equal')
    
    plt.suptitle("Сравнение аналитического и треугольного решений", fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.show()
    
    # Отдельная фигура с подробной статистикой разности
    fig2, ((ax7, ax8), (ax9, ax10)) = plt.subplots(2, 2, figsize=(15, 12))
    
    # 7. Абсолютная разность с контурами
    im7 = ax7.contourf(X1, Y1, Z_diff, levels=50, cmap='RdBu_r')
    contour7 = ax7.contour(X1, Y1, Z_diff, levels=10, colors='black', linewidths=0.5, alpha=0.7)
    ax7.clabel(contour7, inline=True, fontsize=8)
    ax7.set_title("Абсолютная разность с контурами", fontsize=12, fontweight='bold')
    ax7.set_xlabel('X')
    ax7.set_ylabel('Y')
    ax7.set_aspect('equal')
    plt.colorbar(im7, ax=ax7, label='Z1 - Z2')
    
    # 8. Относительная разность с контурами
    im8 = ax8.contourf(X1, Y1, Z_rel_diff_plot, levels=50, cmap='YlOrRd')
    contour8 = ax8.contour(X1, Y1, Z_rel_diff_plot, levels=10, colors='black', linewidths=0.5, alpha=0.7)
    ax8.clabel(contour8, inline=True, fontsize=8)
    ax8.set_title("Относительная разность (%) с контурами", fontsize=12, fontweight='bold')
    ax8.set_xlabel('X')
    ax8.set_ylabel('Y')
    ax8.set_aspect('equal')
    plt.colorbar(im8, ax=ax8, label='Разность (%)')
    
    # 9. Гистограмма абсолютной разности
    diff_flat = Z_diff.flatten()
    ax9.hist(diff_flat, bins=50, edgecolor='black', alpha=0.7, color='blue')
    ax9.axvline(x=0, color='r', linestyle='--', linewidth=2, label='Нулевая разность')
    ax9.set_title("Распределение абсолютной разности", fontsize=12, fontweight='bold')
    ax9.set_xlabel('Z1 - Z2')
    ax9.set_ylabel('Частота')
    ax9.legend()
    ax9.grid(True, alpha=0.3)
    
    # 10. Гистограмма относительной разности (только конечные значения)
    rel_diff_flat = Z_rel_diff.flatten()
    finite_mask = np.isfinite(rel_diff_flat)
    if np.any(finite_mask):
        finite_vals = rel_diff_flat[finite_mask]
        ax10.hist(finite_vals, bins=50, edgecolor='black', alpha=0.7, color='orange')
        ax10.set_title("Распределение относительной разности (%)", fontsize=12, fontweight='bold')
        ax10.set_xlabel('Относительная разность (%)')
        ax10.set_ylabel('Частота')
        ax10.grid(True, alpha=0.3)
    else:
        ax10.text(0.5, 0.5, 'Все значения нулевые', 
                 horizontalalignment='center', verticalalignment='center',
                 transform=ax10.transAxes, fontsize=14)
    
    plt.suptitle("Статистический анализ разности решений", fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.show()
    
    # Третья фигура: 3D поверхности для визуализации разности
    fig3 = plt.figure(figsize=(16, 6))
    
    # 3D поверхность Z1
    ax11 = fig3.add_subplot(1, 3, 1, projection='3d')
    surf1 = ax11.plot_surface(X1, Y1, Z1, cmap='viridis', alpha=0.8, edgecolor='none', 
                              rstride=1, cstride=1, antialiased=True)
    ax11.set_title(f"{title1} - 3D", fontsize=12, fontweight='bold')
    ax11.set_xlabel('X')
    ax11.set_ylabel('Y')
    ax11.set_zlabel('Z')
    fig3.colorbar(surf1, ax=ax11, shrink=0.5, label='Значение')
    
    # 3D поверхность Z2
    ax12 = fig3.add_subplot(1, 3, 2, projection='3d')
    surf2 = ax12.plot_surface(X2, Y2, Z2, cmap='viridis', alpha=0.8, edgecolor='none',
                              rstride=1, cstride=1, antialiased=True)
    ax12.set_title(f"{title2} - 3D", fontsize=12, fontweight='bold')
    ax12.set_xlabel('X')
    ax12.set_ylabel('Y')
    ax12.set_zlabel('Z')
    fig3.colorbar(surf2, ax=ax12, shrink=0.5, label='Значение')
    
    # 3D поверхность разности
    ax13 = fig3.add_subplot(1, 3, 3, projection='3d')
    surf3 = ax13.plot_surface(X1, Y1, Z_diff, cmap='coolwarm', alpha=0.8, edgecolor='none',
                              rstride=1, cstride=1, antialiased=True)
    ax13.set_title("Абсолютная разность (3D)", fontsize=12, fontweight='bold')
    ax13.set_xlabel('X')
    ax13.set_ylabel('Y')
    ax13.set_zlabel('Z1 - Z2')
    fig3.colorbar(surf3, ax=ax13, shrink=0.5, label='Разность')
    
    plt.suptitle("3D визуализация решений и их разности", fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.show()
    
    # Вывод статистики
    print("=" * 70)
    print("СТАТИСТИКА СРАВНЕНИЯ")
    print("=" * 70)
    print(f"Размер сетки: {Z1.shape[0]}x{Z1.shape[1]}")
    print(f"Диапазон координат X: [{ux1[0]:.3e}, {ux1[-1]:.3e}]")
    print(f"Диапазон координат Y: [{uy1[0]:.3e}, {uy1[-1]:.3e}]")
    print(f"Диапазон значений Z1: [{np.min(Z1):.3e}, {np.max(Z1):.3e}]")
    print(f"Диапазон значений Z2: [{np.min(Z2):.3e}, {np.max(Z2):.3e}]")
    print("-" * 70)
    print(f"Максимальная абсолютная разность: {np.max(np.abs(Z_diff)):.6e}")
    print(f"Средняя абсолютная разность: {np.mean(np.abs(Z_diff)):.6e}")
    
    # Для относительной разности считаем только конечные значения
    finite_rel_diff = Z_rel_diff[np.isfinite(Z_rel_diff)]
    if len(finite_rel_diff) > 0:
        print(f"Средняя относительная разность: {np.mean(finite_rel_diff):.2f}%")
        print(f"Максимальная относительная разность: {np.max(finite_rel_diff):.2f}%")
    else:
        print(f"Средняя относительная разность: 0.00%")
        print(f"Максимальная относительная разность: 0.00%")
    
    print(f"RMSE (среднеквадратичная ошибка): {np.sqrt(np.mean(Z_diff**2)):.6e}")
    print(f"Коэффициент корреляции: {np.corrcoef(Z1.flatten(), Z2.flatten())[0,1]:.6f}")
    
    # Проверка на полное совпадение
    if np.allclose(Z1, Z2, atol=1e-12):
        print("-" * 70)
        print("ВНИМАНИЕ: Решения полностью совпадают!")
        print("Абсолютная разность равна нулю для всех точек.")
    else:
        # Находим точки с максимальной разностью
        max_diff_idx = np.unravel_index(np.argmax(np.abs(Z_diff)), Z_diff.shape)
        max_x = X1[max_diff_idx]
        max_y = Y1[max_diff_idx]
        max_z1 = Z1[max_diff_idx]
        max_z2 = Z2[max_diff_idx]
        print("-" * 70)
        print(f"Точка максимальной разности: x={max_x:.3e}, y={max_y:.3e}")
        print(f"Значения: Z1={max_z1:.3e}, Z2={max_z2:.3e}")
    
    print("=" * 70)

# Использование функции
if __name__ == "__main__":
    # Укажите пути к вашим файлам
    file1 = "analytical_solution.txt"
    file2 = "solution_triangular.txt"
    
    try:
        # Построение графиков с разностью
        plot_comparison_with_difference(file1, file2)
    except FileNotFoundError as e:
        print(f"Ошибка: Файл не найден: {e}")
    except Exception as e:
        print(f"Ошибка при построении графиков: {e}")
        import traceback
        traceback.print_exc()