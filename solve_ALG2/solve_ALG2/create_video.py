import cv2
import glob
import re
import os

def extract_number(filename):
    match = re.search(r'v_iteration_(\d+)_contours', filename)
    return int(match.group(1)) if match else 0

files = sorted(glob.glob("v_iteration_*_contours.png"), key=extract_number)
print(f"Найдено файлов: {len(files)}")

if len(files) <= 1:
    print("Недостаточно файлов!")
    # Покажем все файлы в папке для диагностики
    print("\nВсе файлы в папке:")
    all_files = sorted(glob.glob("*.png"))
    for f in all_files:
        print(f"  {f}")
    exit()

# Читаем первый кадр
first_frame = cv2.imread(files[0])
if first_frame is None:
    print("Ошибка чтения первого файла!")
    exit()

height, width = first_frame.shape[:2]
print(f"Размер кадра: {width}x{height}")

# Создаем видео
fps = 2
video_name = 'velocity_animation.mp4'
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
out = cv2.VideoWriter(video_name, fourcc, fps, (width, height))

if not out.isOpened():
    print("Ошибка создания видео!")
    exit()

# Добавляем кадры
for i, file in enumerate(files, 1):
    frame = cv2.imread(file)
    if frame is not None:
        # Изменяем размер если нужно
        if frame.shape[:2] != (height, width):
            frame = cv2.resize(frame, (width, height))
        
        # Добавляем текст с номером итерации
        iteration = extract_number(file)
        cv2.putText(frame, f'Iteration: {iteration}', (10, 30),
                   cv2.FONT_HERSHEY_SIMPLEX, 1, (255, 255, 255), 2)
        
        out.write(frame)
        print(f"Добавлен кадр {i}: {file}")
    else:
        print(f"Ошибка чтения: {file}")

out.release()
print(f"\nВидео создано: {video_name}")
print(f"Всего кадров: {len(files)}")
