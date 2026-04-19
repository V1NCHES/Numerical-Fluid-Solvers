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

def plot_comparison(file1, file2, title1="Аналитическое решение", title2="Треугольное решение"):
    """Строит тепловые карты и контуры для двух файлов"""
    
    # Загружаем данные
    X1, Y1, Z1 = load_data_from_file(file1)
    X2, Y2, Z2 = load_data_from_file(file2)
    
    # Создаем фигуру с 4 подграфиками
    fig = plt.figure(figsize=(16, 12))
    
    # 1. Тепловая карта для файла 1
    ax1 = fig.add_subplot(2, 2, 1)
    contour1 = ax1.contourf(X1, Y1, Z1, levels=50, cmap='viridis')
    ax1.set_title(title1 + " - Тепловая карта")
    ax1.set_xlabel('X')
    ax1.set_ylabel('Y')
    plt.colorbar(contour1, ax=ax1, label='Значение')
    ax1.set_aspect('equal')
    
    # 2. Контуры для файла 1
    ax2 = fig.add_subplot(2, 2, 2)
    contour2 = ax2.contour(X1, Y1, Z1, levels=15, colors='black', linewidths=0.5)
    ax2.clabel(contour2, inline=True, fontsize=8)
    ax2.set_title(title1 + " - Контуры")
    ax2.set_xlabel('X')
    ax2.set_ylabel('Y')
    ax2.set_aspect('equal')
    
    # 3. Тепловая карта для файла 2
    ax3 = fig.add_subplot(2, 2, 3)
    contour3 = ax3.contourf(X2, Y2, Z2, levels=50, cmap='plasma')
    ax3.set_title(title2 + " - Тепловая карта")
    ax3.set_xlabel('X')
    ax3.set_ylabel('Y')
    plt.colorbar(contour3, ax=ax3, label='Значение')
    ax3.set_aspect('equal')
    
    # 4. Контуры для файла 2
    ax4 = fig.add_subplot(2, 2, 4)
    contour4 = ax4.contour(X2, Y2, Z2, levels=15, colors='black', linewidths=0.5)
    ax4.clabel(contour4, inline=True, fontsize=8)
    ax4.set_title(title2 + " - Контуры")
    ax4.set_xlabel('X')
    ax4.set_ylabel('Y')
    ax4.set_aspect('equal')
    
    plt.tight_layout()
    plt.show()
    
    # Также создадим отдельную фигуру для сравнения обоих решений
    fig2, (ax5, ax6) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Сравнение тепловых карт
    im1 = ax5.contourf(X1, Y1, Z1, levels=50, cmap='viridis')
    ax5.set_title(title1)
    ax5.set_xlabel('X')
    ax5.set_ylabel('Y')
    ax5.set_aspect('equal')
    plt.colorbar(im1, ax=ax5, label='Значение')
    
    im2 = ax6.contourf(X2, Y2, Z2, levels=50, cmap='viridis')
    ax6.set_title(title2)
    ax6.set_xlabel('X')
    ax6.set_ylabel('Y')
    ax6.set_aspect('equal')
    plt.colorbar(im2, ax=ax6, label='Значение')
    
    plt.tight_layout()
    plt.show()

# Использование функции
if __name__ == "__main__":
    # Укажите пути к вашим файлам
    file1 = "analytical_solution.txt"
    file2 = "solution_triangular.txt"
    
    # Построение графиков
    plot_comparison(file1, file2)