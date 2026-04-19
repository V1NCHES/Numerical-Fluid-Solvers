import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import re
import argparse
from matplotlib.animation import FuncAnimation

def organize_by_data_type():
    """
    Организует файлы по типам данных
    Возвращает словарь с файлами для каждого типа
    """
    files = glob.glob("*.txt")
    
    organized = {
        'v': [],      # Скорость
        'q1': [],     # Поток q1
        'lambda1': [], # Множитель Лагранжа
        'norm': [],   # Нормы
        'diff': [],   # Разности (diff.txt)
        'other': []   # Остальные файлы
    }
    
    for file in files:
        filename = os.path.basename(file).lower()
        
        if 'v_iteration' in filename:
            organized['v'].append(file)
        elif 'q1_iteration' in filename:
            organized['q1'].append(file)
        elif 'lambda1_iteration' in filename:
            organized['lambda1'].append(file)
        elif 'norm' in filename and 'iteration' in filename:
            organized['norm'].append(file)
        elif 'diff' in filename:
            organized['diff'].append(file)
        else:
            organized['other'].append(file)
    
    return organized

def create_folders():
    """Создает структуру папок для организации графиков"""
    folders = [
        'plots',
        'plots/v_velocity',
        'plots/q1_flow',
        'plots/lambda1_lagrange',
        'plots/norms',
        'plots/v_differences',  # Папка для разностей скорости
        'plots/other',
        'plots/comparisons'
    ]
    
    for folder in folders:
        os.makedirs(folder, exist_ok=True)
    
    return folders

def sort_iteration_files(files):
    """Сортирует файлы по номеру итерации"""
    def get_iteration_number(filename):
        match = re.search(r'iteration_(\d+)', filename)
        return int(match.group(1)) if match else 0
    
    return sorted(files, key=get_iteration_number)

def process_diff_file(file_path):
    """
    Обрабатывает файл diff.txt и строит график сходимости
    """
    try:
        print(f"\nОбработка файла diff.txt...")
        
        # Загружаем данные
        diff_data = np.loadtxt(file_path)
        
        # Номера итераций (начиная с 1)
        iterations = list(range(1, len(diff_data) + 1))
        
        # Создаем несколько графиков
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # 1. График сходимости в log scale
        axes[0].semilogy(iterations, diff_data, 'bo-', linewidth=2, markersize=6)
        axes[0].set_title('Сходимость λ: max|λ_{n+1} - λ_n|')
        axes[0].set_xlabel('Итерация n')
        axes[0].set_ylabel('Максимальная разность (log scale)')
        axes[0].grid(True, alpha=0.3)
        axes[0].set_yscale('log')
        axes[0].axhline(y=1e-6, color='r', linestyle='--', alpha=0.5, label='1e-6')
        axes[0].axhline(y=1e-8, color='g', linestyle='--', alpha=0.5, label='1e-8')
        axes[0].legend()
        
        # 2. График в линейном масштабе
        axes[1].plot(iterations, diff_data, 'ro-', linewidth=2, markersize=6)
        axes[1].set_title('Сходимость λ: max|λ_{n+1} - λ_n|')
        axes[1].set_xlabel('Итерация n')
        axes[1].set_ylabel('Максимальная разность')
        axes[1].grid(True, alpha=0.3)
        
        plt.suptitle('АНАЛИЗ СХОДИМОСТИ МНОЖИТЕЛЕЙ ЛАГРАНЖА (λ)', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        plt.savefig('plots/lambda_convergence.png', 
                   dpi=150, bbox_inches='tight')
        plt.close()
        
        # Простой график для быстрого просмотра
        plt.figure(figsize=(8, 5))
        plt.semilogy(iterations, diff_data, 'bo-', linewidth=2, markersize=8)
        plt.title('Сходимость множителей Лагранжа')
        plt.xlabel('Итерация')
        plt.ylabel('max|λ_{n+1} - λ_n| (log scale)')
        plt.grid(True, alpha=0.3)
        plt.axhline(y=1e-6, color='r', linestyle='--', alpha=0.7, label='Порог 1e-6')
        plt.legend()
        plt.tight_layout()
        plt.savefig('plots/lambda_convergence_simple.png', 
                   dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"  Графики сходимости сохранены в plots/")
        
        return True
        
    except Exception as e:
        print(f"  Ошибка при обработке diff.txt: {e}")
        return False

def visualize_velocity(file_path, output_folder, x_min=0, x_max=1, y_min=0, y_max=1):
    """
    Визуализирует файл скорости (v) - только 3D и контуры
    """
    try:
        # Загружаем данные
        data = np.loadtxt(file_path)
        filename = os.path.basename(file_path)
        name = os.path.splitext(filename)[0]
        
        print(f"  Визуализация скорости: {filename} (размер: {data.shape})")
        
        # Определяем размеры
        ny, nx = data.shape
        
        # Создаем координатную сетку
        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_min, y_max, ny)
        X, Y = np.meshgrid(x, y)
        
        # 1. 3D ПОВЕРХНОСТЬ
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Для больших сеток используем downsampling для 3D графика
        stride_x = max(1, nx // 50)
        stride_y = max(1, ny // 50)
        
        surf = ax.plot_surface(X[::stride_y, ::stride_x], Y[::stride_y, ::stride_x], 
                              data[::stride_y, ::stride_x],
                              cmap='viridis', alpha=0.8,
                              rstride=1, cstride=1, linewidth=0)
        
        ax.set_title(f'{name}\nСкорость (v) - 3D')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Скорость')
        ax.view_init(elev=30, azim=45)
        
        plt.colorbar(surf, shrink=0.6, pad=0.1, label='Скорость')
        plt.tight_layout()
        plt.savefig(f'{output_folder}/{name}_3d.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        # 2. КОНТУРНЫЙ ГРАФИК
        plt.figure(figsize=(10, 8))
        
        # Контурная диаграмма
        contour = plt.contourf(X, Y, data, levels=30, cmap='viridis')
        plt.contour(X, Y, data, levels=15, colors='black', linewidths=0.5, alpha=0.5)
        
        plt.title(f'{name}\nСкорость (v) - Контуры')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.colorbar(contour, label='Скорость')
        plt.tight_layout()
        plt.savefig(f'{output_folder}/{name}_contours.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        return True
        
    except Exception as e:
        print(f"  Ошибка при обработке {file_path}: {e}")
        return False

def visualize_other(file_path, output_folder, x_min=0, x_max=1, y_min=0, y_max=1):
    """
    Визуализирует другие файлы (q1, lambda1, norm) - только контуры и тепловую карту
    """
    try:
        # Загружаем данные
        data = np.loadtxt(file_path)
        filename = os.path.basename(file_path)
        name = os.path.splitext(filename)[0]
        
        print(f"  Визуализация: {filename} (размер: {data.shape})")
        
        # Определяем тип данных из имени файла
        if 'q1' in filename.lower():
            data_type = 'Поток (q1)'
            cmap = 'coolwarm'
        elif 'lambda' in filename.lower():
            data_type = 'Множитель Лагранжа (λ)'
            cmap = 'RdBu_r'
        elif 'norm' in filename.lower():
            data_type = 'Норма'
            cmap = 'plasma'
        else:
            data_type = 'Данные'
            cmap = 'viridis'
        
        # Определяем размеры
        ny, nx = data.shape
        
        # Создаем координатную сетку
        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_min, y_max, ny)
        X, Y = np.meshgrid(x, y)
        
        # 1. ТЕПЛОВАЯ КАРТА
        plt.figure(figsize=(10, 8))
        im = plt.imshow(data, cmap=cmap, origin='lower',
                       extent=[x_min, x_max, y_min, y_max], aspect='auto')
        
        plt.title(f'{name}\n{data_type} - Тепловая карта')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.colorbar(im, label=data_type.split('(')[-1].replace(')', '') if '(' in data_type else 'Значение')
        plt.tight_layout()
        plt.savefig(f'{output_folder}/{name}_heatmap.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        # 2. КОНТУРНЫЙ ГРАФИК
        plt.figure(figsize=(10, 8))
        
        # Контурная диаграмма
        contour = plt.contourf(X, Y, data, levels=30, cmap=cmap)
        plt.contour(X, Y, data, levels=15, colors='black', linewidths=0.5, alpha=0.5)
        
        plt.title(f'{name}\n{data_type} - Контуры')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.colorbar(contour, label=data_type.split('(')[-1].replace(')', '') if '(' in data_type else 'Значение')
        plt.tight_layout()
        plt.savefig(f'{output_folder}/{name}_contours.png', dpi=150, bbox_inches='tight')
        plt.close()
        
        return True
        
    except Exception as e:
        print(f"  Ошибка при обработке {file_path}: {e}")
        return False

def calculate_velocity_differences(v_files, x_min=0, x_max=1, y_min=0, y_max=1):
    """
    Вычисляет и визуализирует разности между последовательными итерациями скорости
    """
    if len(v_files) < 2:
        print("\nНедостаточно файлов скорости для вычисления разностей")
        return
    
    print(f"\n{'='*60}")
    print("ВЫЧИСЛЕНИЕ РАЗНОСТЕЙ МЕЖДУ ИТЕРАЦИЯМИ СКОРОСТИ")
    print(f"{'='*60}")
    
    # Сортируем файлы по итерациям
    sorted_files = sort_iteration_files(v_files)
    
    # Создаем папку для разностей
    diff_folder = 'plots/v_differences'
    os.makedirs(diff_folder, exist_ok=True)
    
    # Загружаем все данные скорости
    velocity_data = []
    iterations = []
    
    for file_path in sorted_files:
        try:
            data = np.loadtxt(file_path)
            velocity_data.append(data)
            
            # Извлекаем номер итерации
            match = re.search(r'iteration_(\d+)', file_path)
            if match:
                iterations.append(int(match.group(1)))
            else:
                iterations.append(len(iterations))
        except Exception as e:
            print(f"  Ошибка загрузки {file_path}: {e}")
            continue
    
    if len(velocity_data) < 2:
        print("  Недостаточно данных для вычисления разностей")
        return
    
    print(f"  Загружено {len(velocity_data)} итераций скорости")
    
    # Создаем список максимальных разностей для графика сходимости
    max_differences = []
    
    # Вычисляем разности между последовательными итерациями
    for i in range(len(velocity_data) - 1):
        try:
            # Вычисляем разность (v_{n+1} - v_n)
            diff = velocity_data[i + 1] - velocity_data[i]
            
            # Максимальная разность по модулю
            max_diff = np.abs(diff).max()
            max_differences.append(max_diff)
            
            # Визуализируем разность (контуры и тепловая карта)
            visualize_difference(diff, iterations[i], iterations[i+1], 
                               diff_folder, x_min, x_max, y_min, y_max)
            
            print(f"  Итерации {iterations[i]} → {iterations[i+1]}: max|Δv| = {max_diff:.6e}")
            
        except Exception as e:
            print(f"  Ошибка вычисления разности {i}→{i+1}: {e}")
    
    # График сходимости разностей скорости
    if max_differences:
        plot_velocity_convergence(max_differences, iterations, diff_folder)
        
        # Сохраняем max_differences в файл
        np.savetxt(os.path.join(diff_folder, 'velocity_max_differences.txt'), 
                  max_differences)
        
        print(f"\n  Максимальные разности сохранены в {diff_folder}/velocity_max_differences.txt")

def visualize_difference(diff_data, iter_from, iter_to, output_folder, 
                        x_min, x_max, y_min, y_max):
    """
    Визуализирует разность между двумя итерациями скорости
    """
    try:
        # Определяем размеры
        ny, nx = diff_data.shape
        
        # Создаем координатную сетку
        x = np.linspace(x_min, x_max, nx)
        y = np.linspace(y_min, y_max, ny)
        X, Y = np.meshgrid(x, y)
        
        # 1. ТЕПЛОВАЯ КАРТА РАЗНОСТИ
        plt.figure(figsize=(10, 8))
        
        # Используем diverging colormap для разностей
        im = plt.imshow(diff_data, cmap='RdBu_r', origin='lower',
                       extent=[x_min, x_max, y_min, y_max], aspect='auto')
        
        plt.title(f'Разность скорости: v_{iter_to} - v_{iter_from}')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.colorbar(im, label='Δv')
        plt.tight_layout()
        
        filename = f'velocity_diff_{iter_from}_to_{iter_to}_heatmap.png'
        plt.savefig(os.path.join(output_folder, filename), 
                   dpi=150, bbox_inches='tight')
        plt.close()
        
        # 2. КОНТУРНЫЙ ГРАФИК РАЗНОСТИ
        plt.figure(figsize=(10, 8))
        
        # Контурная диаграмма
        contour = plt.contourf(X, Y, diff_data, levels=30, cmap='RdBu_r')
        plt.contour(X, Y, diff_data, levels=15, colors='black', 
                   linewidths=0.5, alpha=0.5)
        
        plt.title(f'Разность скорости: v_{iter_to} - v_{iter_from}')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.colorbar(contour, label='Δv')
        plt.tight_layout()
        
        filename = f'velocity_diff_{iter_from}_to_{iter_to}_contours.png'
        plt.savefig(os.path.join(output_folder, filename), 
                   dpi=150, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"    Ошибка визуализации разности {iter_from}→{iter_to}: {e}")

def plot_velocity_convergence(max_differences, iterations, output_folder):
    """
    Строит график сходимости разностей скорости
    """
    try:
        # Итерации для разностей (между итерациями)
        diff_iterations = iterations[1:]  # Начиная со второй итерации
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # 1. График сходимости в log scale
        axes[0].semilogy(diff_iterations, max_differences, 'bo-', 
                       linewidth=2, markersize=6)
        axes[0].set_title('Сходимость скорости: max|v_{n+1} - v_n|')
        axes[0].set_xlabel('Итерация n')
        axes[0].set_ylabel('Максимальная разность (log scale)')
        axes[0].grid(True, alpha=0.3)
        axes[0].axhline(y=1e-6, color='r', linestyle='--', alpha=0.5, label='1e-6')
        axes[0].axhline(y=1e-8, color='g', linestyle='--', alpha=0.5, label='1e-8')
        axes[0].legend()
        
        # 2. График в линейном масштабе
        axes[1].plot(diff_iterations, max_differences, 'ro-', 
                    linewidth=2, markersize=6)
        axes[1].set_title('Сходимость скорости: max|v_{n+1} - v_n|')
        axes[1].set_xlabel('Итерация n')
        axes[1].set_ylabel('Максимальная разность')
        axes[1].grid(True, alpha=0.3)
        
        plt.suptitle('СХОДИМОСТЬ СКОРОСТИ (v)', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        plt.savefig(os.path.join(output_folder, 'velocity_convergence.png'), 
                   dpi=150, bbox_inches='tight')
        plt.close()
        
        # Простой график для быстрого просмотра
        plt.figure(figsize=(8, 5))
        plt.semilogy(diff_iterations, max_differences, 'bo-', 
                    linewidth=2, markersize=8)
        plt.title('Сходимость скорости: max|v_{n+1} - v_n|')
        plt.xlabel('Итерация')
        plt.ylabel('Максимальная разность (log scale)')
        plt.grid(True, alpha=0.3)
        plt.axhline(y=1e-6, color='r', linestyle='--', alpha=0.7, label='Порог 1e-6')
        plt.legend()
        plt.tight_layout()
        
        plt.savefig(os.path.join(output_folder, 'velocity_convergence_simple.png'), 
                   dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"  Графики сходимости скорости сохранены в {output_folder}/")
        
    except Exception as e:
        print(f"  Ошибка построения графика сходимости скорости: {e}")

def visualize_by_type(data_type, files, x_min=0, x_max=1, y_min=0, y_max=1):
    """
    Визуализирует все файлы определенного типа
    """
    if not files:
        print(f"Нет файлов для типа: {data_type}")
        return 0
    
    # Определяем папку для этого типа данных
    type_folders = {
        'v': 'plots/v_velocity',
        'q1': 'plots/q1_flow',
        'lambda1': 'plots/lambda1_lagrange',
        'norm': 'plots/norms',
        'diff': 'plots',  # diff идет в корень plots
        'other': 'plots/other'
    }
    
    output_folder = type_folders.get(data_type, 'plots/other')
    os.makedirs(output_folder, exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"ТИП ДАННЫХ: {data_type.upper()}")
    print(f"Количество файлов: {len(files)}")
    print(f"Папка для графиков: {output_folder}")
    print(f"{'='*60}")
    
    # Особый случай для diff файлов
    if data_type == 'diff':
        success_count = 0
        for file_path in files:
            if process_diff_file(file_path):
                success_count += 1
        return success_count
    
    # Сортируем файлы по итерациям
    sorted_files = sort_iteration_files(files)
    
    success_count = 0
    for file_path in sorted_files:
        if data_type == 'v':
            # Для скорости - только 3D и контуры
            if visualize_velocity(file_path, output_folder, x_min, x_max, y_min, y_max):
                success_count += 1
        else:
            # Для остальных - только тепловая карта и контуры
            if visualize_other(file_path, output_folder, x_min, x_max, y_min, y_max):
                success_count += 1
    
    print(f"\nУспешно обработано: {success_count}/{len(files)} файлов")
    
    # Для скорости дополнительно вычисляем разности
    if data_type == 'v' and len(sorted_files) > 1:
        calculate_velocity_differences(sorted_files, x_min, x_max, y_min, y_max)
    
    return success_count

def create_summary_report(organized_files, total_processed):
    """
    Создает сводный отчет
    """
    print(f"\n{'='*60}")
    print("СВОДНЫЙ ОТЧЕТ")
    print(f"{'='*60}")
    
    for data_type, files in organized_files.items():
        if files and data_type != 'other':
            type_name = {
                'v': 'Скорость (v)',
                'q1': 'Поток (q1)',
                'lambda1': 'Множители Лагранжа',
                'norm': 'Нормы',
                'diff': 'Разности λ'
            }.get(data_type, data_type)
            
            print(f"{type_name:25} : {len(files):3d} файлов")
    
    print(f"{'='*60}")
    print(f"Всего обработано файлов: {total_processed}")

def main():
    parser = argparse.ArgumentParser(description='Визуализация данных (упрощенная версия)')
    parser.add_argument('--xmin', type=float, default=0.0, help='Минимальное значение X')
    parser.add_argument('--xmax', type=float, default=1.0, help='Максимальное значение X')
    parser.add_argument('--ymin', type=float, default=0.0, help='Минимальное значение Y')
    parser.add_argument('--ymax', type=float, default=1.0, help='Максимальное значение Y')
    parser.add_argument('--type', type=str, 
                       choices=['all', 'v', 'q1', 'lambda1', 'norm', 'diff'],
                       default='all', help='Тип данных для визуализации')
    
    args = parser.parse_args()
    
    print("="*70)
    print("ВИЗУАЛИЗАЦИЯ ДАННЫХ (УПРОЩЕННАЯ ВЕРСИЯ)")
    print("="*70)
    print(f"Расчетная область: X=[{args.xmin:.2f}, {args.xmax:.2f}], Y=[{args.ymin:.2f}, {args.ymax:.2f}]")
    print(f"Тип данных: {args.type}")
    print("="*70)
    print("РЕЖИМЫ ВИЗУАЛИЗАЦИИ:")
    print("  • Скорость (v): 3D поверхность + контуры")
    print("  • Поток (q1): тепловая карта + контуры")
    print("  • Множители Лагранжа: тепловая карта + контуры")
    print("  • Нормы: тепловая карта + контуры")
    print("  • Разности скорости: тепловая карта + контуры")
    print("  • Разности λ: графики сходимости")
    print("="*70)
    
    # Создаем структуру папок
    create_folders()
    
    # Организуем файлы по типам
    organized_files = organize_by_data_type()
    
    # Визуализируем выбранные типы данных
    total_processed = 0
    
    if args.type == 'all':
        # Обрабатываем все типы в определенном порядке
        processing_order = ['diff', 'v', 'q1', 'lambda1', 'norm']
        
        for data_type in processing_order:
            count = visualize_by_type(data_type, organized_files[data_type], 
                                    args.xmin, args.xmax, args.ymin, args.ymax)
            total_processed += count
    else:
        # Обрабатываем только выбранный тип
        count = visualize_by_type(args.type, organized_files[args.type],
                                args.xmin, args.xmax, args.ymin, args.ymax)
        total_processed = count
    
    # Создаем сводный отчет
    create_summary_report(organized_files, total_processed)
    
    print(f"\n{'='*70}")
    print("🎉 ВИЗУАЛИЗАЦИЯ ЗАВЕРШЕНА!")
    print(f"{'='*70}")
    print("\n📁 СТРУКТУРА ПАПОК:")
    print("  plots/v_velocity/       - Скорость (3D + контуры)")
    print("  plots/q1_flow/          - Поток (тепловая + контуры)")
    print("  plots/lambda1_lagrange/ - Множители Лагранжа (тепловая + контуры)")
    print("  plots/norms/            - Нормы (тепловая + контуры)")
    print("  plots/v_differences/    - Разности скорости (тепловая + контуры)")
    print("\n📊 КЛЮЧЕВЫЕ ГРАФИКИ:")
    print("  plots/lambda_convergence.png - Сходимость λ")
    print("  plots/v_differences/velocity_convergence.png - Сходимость скорости")
    print(f"\n📈 ОБРАБОТАНО ФАЙЛОВ: {total_processed}")
    print(f"{'='*70}")

if __name__ == "__main__":
    main()