import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import re
import argparse
from matplotlib.animation import FuncAnimation
from PIL import Image
import imageio

# [Все существующие функции остаются без изменений до конца файла]
# ... (весь ваш существующий код до main)

def create_velocity_animation(v_files, grid_size_x=1.0, grid_size_y=1.0, 
                             output_gif='velocity_animation.gif', fps=5, 
                             vmin=None, vmax=None, cmap='viridis'):
    """
    Создает анимацию из файлов скорости
    
    Параметры:
    - v_files: список файлов скорости
    - grid_size_x, grid_size_y: размеры области
    - output_gif: имя выходного GIF файла
    - fps: кадров в секунду
    - vmin, vmax: границы цветовой шкалы (если None, определяются автоматически)
    - cmap: цветовая карта
    """
    if not v_files:
        print("Нет файлов скорости для создания анимации")
        return
    
    print(f"\n{'='*60}")
    print("СОЗДАНИЕ АНИМАЦИИ СКОРОСТИ")
    print(f"{'='*60}")
    
    # Сортируем файлы по итерациям
    sorted_files = sort_iteration_files(v_files)
    print(f"Найдено файлов: {len(sorted_files)}")
    
    # Загружаем все данные для определения глобальных границ
    all_data = []
    data_list = []
    x_coords_list = []
    y_coords_list = []
    iterations = []
    
    for file_path in sorted_files:
        data, x_coords, y_coords = load_data_from_txt(file_path, grid_size_x, grid_size_y)
        if data is not None:
            data_list.append(data)
            all_data.append(data.flatten())
            x_coords_list.append(x_coords)
            y_coords_list.append(y_coords)
            
            # Извлекаем номер итерации
            match = re.search(r'iteration_(\d+)', file_path)
            iterations.append(int(match.group(1)) if match else len(iterations))
    
    if not data_list:
        print("Не удалось загрузить данные")
        return
    
    # Определяем границы цветовой шкалы
    if vmin is None:
        vmin = np.min(all_data)
    if vmax is None:
        vmax = np.max(all_data)
    
    print(f"Границы цветовой шкалы: vmin={vmin:.6f}, vmax={vmax:.6f}")
    
    # Создаем фигуру для анимации
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Берем координаты из первого файла
    X, Y = np.meshgrid(x_coords_list[0], y_coords_list[0])
    
    # Создаем график с фиксированной цветовой шкалой
    im = ax.pcolormesh(X, Y, data_list[0], vmin=vmin, vmax=vmax, 
                      cmap=cmap, shading='auto')
    
    # Добавляем colorbar
    cbar = plt.colorbar(im, ax=ax, label='Скорость (v)')
    
    # Настраиваем оси
    ax.set_xlabel('X координата')
    ax.set_ylabel('Y координата')
    ax.set_xlim(0, grid_size_x)
    ax.set_ylim(0, grid_size_y)
    
    # Заголовок с номером итерации
    title = ax.set_title(f'Скорость - Итерация: {iterations[0]}')
    
    # Функция обновления для анимации
    def update(frame):
        # Обновляем данные
        im.set_array(data_list[frame].ravel())
        
        # Обновляем заголовок
        title.set_text(f'Скорость - Итерация: {iterations[frame]}')
        
        return [im, title]
    
    # Создаем анимацию
    anim = FuncAnimation(
        fig, update, frames=len(data_list),
        interval=1000/fps,  # интервал в миллисекундах
        blit=True,
        repeat=True
    )
    
    # Сохраняем анимацию
    output_path = os.path.join('plots', output_gif)
    print(f"Сохранение анимации в {output_path}...")
    
    # Пробуем сохранить как GIF
    try:
        anim.save(output_path, writer='pillow', fps=fps, dpi=100)
        print(f"✅ Анимация успешно сохранена: {output_path}")
    except Exception as e:
        print(f"❌ Ошибка при сохранении GIF: {e}")
        # Пробуем альтернативный метод
        try:
            # Сохраняем как MP4
            mp4_path = output_path.replace('.gif', '.mp4')
            anim.save(mp4_path, writer='ffmpeg', fps=fps, dpi=100)
            print(f"✅ Анимация сохранена как MP4: {mp4_path}")
        except:
            print("❌ Не удалось сохранить анимацию")
    
    plt.close(fig)
    
    # Также создаем упрощенную анимацию из изображений (если уже есть сохраненные графики)
    create_animation_from_images(iterations, vmin, vmax)

def create_animation_from_images(iterations, vmin, vmax):
    """
    Создает анимацию из уже сохраненных изображений скорости
    """
    # Ищем все 3D графики скорости
    image_files = glob.glob('plots/v_velocity/*iteration_*_3d.png')
    
    if not image_files:
        print("Не найдены сохраненные изображения для анимации")
        return
    
    # Сортируем изображения по итерациям
    def get_iter_num(filename):
        match = re.search(r'iteration_(\d+)', filename)
        return int(match.group(1)) if match else 0
    
    image_files = sorted(image_files, key=get_iter_num)
    
    if len(image_files) < 2:
        print("Недостаточно изображений для анимации")
        return
    
    print(f"\nСоздание анимации из {len(image_files)} изображений...")
    
    # Создаем GIF из изображений
    gif_path = 'plots/velocity_animation_from_images.gif'
    
    try:
        with imageio.get_writer(gif_path, mode='I', fps=2) as writer:
            for img_file in image_files:
                image = imageio.imread(img_file)
                writer.append_data(image)
        
        print(f"✅ Анимация из изображений сохранена: {gif_path}")
    except Exception as e:
        print(f"❌ Ошибка при создании анимации из изображений: {e}")

def add_animation_to_main():
    """
    Функция для добавления параметров анимации в основной скрипт
    """
    # Добавляем новые аргументы командной строки
    parser = argparse.ArgumentParser(description='Визуализация данных с анимацией')
    parser.add_argument('--sizex', type=float, default=2.0, 
                       help='Размер области по оси X')
    parser.add_argument('--sizey', type=float, default=2.0, 
                       help='Размер области по оси Y')
    parser.add_argument('--type', type=str, 
                       choices=['all', 'v', 'q1', 'lambda1', 'norm', 'norm_lambda1', 'diff'],
                       default='all', help='Тип данных для визуализации')
    
    # НОВЫЕ АРГУМЕНТЫ ДЛЯ АНИМАЦИИ
    parser.add_argument('--animate', action='store_true', 
                       help='Создать анимацию скорости')
    parser.add_argument('--fps', type=int, default=5,
                       help='Количество кадров в секунду для анимации')
    parser.add_argument('--vmin', type=float, default=None,
                       help='Минимальное значение для цветовой шкалы')
    parser.add_argument('--vmax', type=float, default=None,
                       help='Максимальное значение для цветовой шкалы')
    parser.add_argument('--cmap', type=str, default='viridis',
                       help='Цветовая карта для анимации')
    
    return parser

# Модифицируем функцию main() для поддержки анимации
def main_with_animation():
    # Используем новый парсер с аргументами анимации
    parser = add_animation_to_main()
    args = parser.parse_args()
    
    print("="*70)
    print("ВИЗУАЛИЗАЦИЯ ДАННЫХ С АНИМАЦИЕЙ")
    print("="*70)
    print(f"Размер области: X=[0, {args.sizex:.2f}], Y=[0, {args.sizey:.2f}]")
    print(f"Тип данных: {args.type}")
    if args.animate:
        print(f"🎬 АНИМАЦИЯ ВКЛЮЧЕНА: fps={args.fps}, cmap={args.cmap}")
    print("="*70)
    
    # Создаем структуру папок
    create_folders()
    
    # Организуем файлы по типам
    organized_files = organize_by_data_type()
    
    # Выводим информацию о найденных файлах
    for data_type, files in organized_files.items():
        if files and data_type != 'other':
            print(f"{data_type:15}: {len(files)} файлов")
    
    # Визуализируем выбранные типы данных
    total_processed = 0
    
    if args.type == 'all':
        # Обрабатываем все типы в определенном порядке
        processing_order = ['diff', 'v', 'q1', 'lambda1', 'norm', 'norm_lambda1']
        
        for data_type in processing_order:
            count = visualize_by_type(data_type, organized_files[data_type], 
                                    args.sizex, args.sizey)
            total_processed += count
    else:
        # Обрабатываем только выбранный тип
        count = visualize_by_type(args.type, organized_files[args.type],
                                args.sizex, args.sizey)
        total_processed = count
    
    # СОЗДАЕМ АНИМАЦИЮ ЕСЛИ ЗАПРОШЕНО
    if args.animate and organized_files['v']:
        print("\n" + "="*60)
        print("ШАГ 2: СОЗДАНИЕ АНИМАЦИИ СКОРОСТИ")
        print("="*60)
        
        create_velocity_animation(
            organized_files['v'],
            grid_size_x=args.sizex,
            grid_size_y=args.sizey,
            output_gif='velocity_animation.gif',
            fps=args.fps,
            vmin=args.vmin,
            vmax=args.vmax,
            cmap=args.cmap
        )
    
    # Создаем сводный отчет
    create_summary_report(organized_files, total_processed)
    
    print(f"\n{'='*70}")
    print("🎉 ВИЗУАЛИЗАЦИЯ ЗАВЕРШЕНА!")
    print(f"{'='*70}")
    
    if args.animate:
        print("\n🎬 АНИМАЦИЯ:")
        print("  - plots/velocity_animation.gif")
        print("  - plots/velocity_animation_from_images.gif (если есть изображения)")
    
    print(f"\n📈 ОБРАБОТАНО ФАЙЛОВ: {total_processed}")
    print(f"{'='*70}")

# Для обратной совместимости оставляем старую main
if __name__ == "__main__":
    # Проверяем, есть ли аргумент --animate
    if '--animate' in sys.argv:
        main_with_animation()
    else:
        main()