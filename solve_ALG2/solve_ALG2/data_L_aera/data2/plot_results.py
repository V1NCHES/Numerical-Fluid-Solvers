import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import argparse
import re
import gc

# Конфигурация по умолчанию
DEFAULT_FOLDER = "output_data"
DEFAULT_OUT_FOLDER = "plots"
DEFAULT_TAU_S_RANGE = [0.0, 0.2, 0.202] # [min, max, threshold]

def load_data(filename):
    try:
        data = np.loadtxt(filename)
        if data.size == 0: return None
        return data
    except Exception as e:
        print(f"Ошибка загрузки {filename}: {e}")
        return None

def plot_convergence(filename, output_folder):
    data = load_data(filename)
    if data is None: return
    
    plt.figure(figsize=(10, 6))
    plt.plot(data, linewidth=2, color='#1f77b4')
    plt.yscale('log')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.xlabel('Iteration')
    plt.ylabel('Max Difference')
    plt.title('Solver Convergence')
    
    os.makedirs(output_folder, exist_ok=True)
    out_path = os.path.join(output_folder, "convergence_log.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Convergence plot saved: {out_path}")

def get_best_colormap(z):
    if z.min() < 0 and z.max() > 0: return 'RdBu_r'
    return 'viridis'

def plot_field(filename, title, output_folder, special_mode=False, 
               show_streamlines=False, show_arcs=False, overlay_rigid=False,
               show_contours=False, tau_range=DEFAULT_TAU_S_RANGE):
    data = load_data(filename)
    if data is None: return

    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    ux, uy = np.sort(np.unique(x)), np.sort(np.unique(y))
    nx, ny = len(ux), len(uy)

    if nx < 2 or ny < 2 or len(z) < 3:
        print(f"Пропуск {filename}: недостаточно точек ({nx}x{ny})")
        return

    # Настройка стиля
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_facecolor('white')
    fig.patch.set_facecolor('white')
    
    plt.rcParams.update({'font.size': 22})
    ax.tick_params(axis='both', which='major', labelsize=24)
    for spine in ax.spines.values():
        spine.set_linewidth(2)

    if not special_mode:
        try: plt.style.use('seaborn-v10_8-darkgrid')
        except: ax.grid(True, linestyle='--', alpha=0.4)
        cmap = get_best_colormap(z)
        levels = 50

    # --- Геометрия L-области для маскирования ---
    ux_unique = np.sort(np.unique(x))
    uy_unique = np.sort(np.unique(y))
    h_x = ux_unique[1] - ux_unique[0] if len(ux_unique) > 1 else 0.01
    h_y = uy_unique[1] - uy_unique[0] if len(uy_unique) > 1 else 0.01
    Lx_total = ux_unique.max() + h_x/2
    Ly_total = uy_unique.max() + h_y/2
    
    bx_detected = x[y > Ly_total * 0.9].max() + h_x/2 if any(y > Ly_total * 0.9) else Lx_total
    by_detected = y[x > Lx_total * 0.9].max() + h_y/2 if any(x > Lx_total * 0.9) else Ly_total
    
    # --- Устранение артефактов на стыке блоков (Point Reconciliation) ---
    z_healed = z.copy()
    # Находим точки с идентичными координатами (граница блоков) и усредняем их значения
    coords = np.column_stack((x, y))
    rounded_coords = np.round(coords, decimals=8)
    unique_coords, inverse_indices, counts = np.unique(rounded_coords, axis=0, return_inverse=True, return_counts=True)
    
    duplicate_indices = np.where(counts > 1)[0]
    for idx_val in duplicate_indices:
        mask_dup = (inverse_indices == idx_val)
        z_healed[mask_dup] = np.mean(z[mask_dup])

    # Создаем триангуляцию и маскируем пустой квадрант
    import matplotlib.tri as tri
    triang = tri.Triangulation(x, y)
    x_mid = x[triang.triangles].mean(axis=1)
    y_mid = y[triang.triangles].mean(axis=1)
    mask = (x_mid > bx_detected + 1e-4) & (y_mid > by_detected + 1e-4)
    triang.set_mask(mask)

    cnt = None
    if special_mode:
        # Заливка жестких зон
        ax.tricontourf(triang, z_healed, levels=[0, tau_range[1]], colors=['#706b67'])
        # Отрисовка границ жестких зон (yield surfaces)
        ax.tricontour(triang, z_healed, levels=[tau_range[1] * 1.001], colors='black', linewidths=1.5)
        
        # --- Скрытие линии на внутреннем стыке блоков (y = by_detected, x < bx) ---
        # Мы перекрываем только внутренний интерфейс белой линией, если там нет жесткой зоны
        ax.plot([0, bx_detected], [by_detected, by_detected], color='white', linewidth=3, zorder=18)
    else:
        limit = max(abs(z.min()), abs(z.max())) or 1.0
        cnt = ax.tricontourf(triang, z_healed, levels=levels, cmap=cmap, vmin=0, vmax=limit)
        if show_contours:
            ax.tricontour(triang, z, levels=15, colors='black', linewidths=0.5, alpha=0.2)

    # Наложение линий тока
    if show_streamlines:
        dir_path = os.path.dirname(filename)
        match = re.search(r'(iter|iret|iteration)_?(\d+)', os.path.basename(filename))
        if match:
            prefix, num = match.group(1), match.group(2)
            psi_filename = None
            for pattern in [f"psi_iteration_{num}.txt", f"{prefix}_{num}_psi.txt", f"psi_{num}.txt"]:
                p_path = os.path.join(dir_path, pattern)
                if os.path.exists(p_path):
                    psi_filename = p_path
                    break
            if psi_filename:
                psi_data = load_data(psi_filename)
                if psi_data is not None:
                    px, py, pz = psi_data[:,0], psi_data[:,1], psi_data[:,2]
                    pux, puy = np.sort(np.unique(px)), np.sort(np.unique(py))
                    if len(pz) == len(pux)*len(puy):
                        ax.contour(pux, puy, pz.reshape((len(puy), len(pux))), 
                                   levels=18, colors='grey', linewidths=1.0, alpha=0.8)

    # Отрисовка теоретических дуг и границ L-области
    if show_arcs or special_mode:
        # --- Отрисовка границы L-области (Черная линия) ---
        # Рисуем только внешние границы, исключая линию y = by_detected
        segments = [
            [(0, 0), (2, 0)],                 # Низ (от 0 до 2)
            [(2, 0), (2, 1)],                 # Правый край (от y=0 до y=1)
            [(1, 1), (2, 1)],                 # Внутр. горизонталь (от x=1 до x=2)
            [(1, 1), (1, 2)],                 # Внутр. вертикаль (от y=1 до y=2)
            [(0, 2), (1, 2)],                 # Верх (от x=0 до x=1)
            [(0, 2), (0, 0)]                  # Лево (от y=0 до y=2)
        ]
        for seg in segments:
            sx, sy = zip(*seg)
            ax.plot(sx, sy, color='black', linewidth=2.5, zorder=20)

        if show_arcs:
            K_const = 1.0 - np.pi / 4.0
            Area_total = Lx_total * by_detected + (Ly_total - by_detected) * bx_detected
            L_val = max(Lx_total, Ly_total)
            discriminant = 16 * L_val**2 - 20 * Area_total * K_const
            l_star = (4 * L_val - np.sqrt(max(0, discriminant))) / (10 * K_const)
            
            R = l_star
            t = np.linspace(0, np.pi/2, 100)
            # 5 вогнутых угловых дуг
            ax.plot(R + R * np.cos(t + np.pi), R + R * np.sin(t + np.pi), color='red', linewidth=4, zorder=25)
            ax.plot(R + R * np.cos(t + np.pi/2), Ly_total - R + R * np.sin(t + np.pi/2), color='red', linewidth=4, zorder=25)
            ax.plot(bx_detected - R + R * np.cos(t), Ly_total - R + R * np.sin(t), color='red', linewidth=4, zorder=25)
            ax.plot(Lx_total - R + R * np.cos(t + 3*np.pi/2), R + R * np.sin(t + 3*np.pi/2), color='red', linewidth=4, zorder=25)
            ax.plot(Lx_total - R + R * np.cos(t), by_detected - R + R * np.sin(t), color='red', linewidth=4, zorder=25)

    if cnt and not special_mode:
        fig.colorbar(cnt).ax.set_ylabel('Value', rotation=270, labelpad=15)
    
    ax.set_aspect('equal')
    plt.tight_layout()

    os.makedirs(output_folder, exist_ok=True)
    suffix = "_rigid_arcs" if show_arcs else ("_yield_sci" if special_mode else "_full")
    out_name = os.path.join(output_folder, os.path.basename(filename).replace('.txt', f"{suffix}.png"))
    
    plt.savefig(out_name, dpi=300)
    plt.close(fig)
    print(f"Saved: {out_name}")
    gc.collect()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=str, default=DEFAULT_FOLDER)
    parser.add_argument('--iter', type=int)
    parser.add_argument('--out', type=str, default=DEFAULT_OUT_FOLDER)
    parser.add_argument('--tau_min', type=float, default=DEFAULT_TAU_S_RANGE[0])
    parser.add_argument('--tau_max', type=float, default=DEFAULT_TAU_S_RANGE[1])
    args = parser.parse_args()
    
    current_tau_range = [args.tau_min, args.tau_max, DEFAULT_TAU_S_RANGE[1]]
    data_folder = args.folder
    if not os.path.exists(data_folder):
        data_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), args.folder)

    if os.path.exists(data_folder):
        # 1. График сходимости
        for f_name in ['diff_max.txt', 'diff_max']:
            p = os.path.join(data_folder, f_name)
            if os.path.exists(p):
                plot_convergence(p, args.out)
                break

        # 2. Список файлов для обработки
        files = sorted([f for f in os.listdir(data_folder) 
                        if f.endswith('.txt') and 'diff_max' not in f])
        
        if args.iter is not None:
            files = [f for f in files if any(int(n) == args.iter for n in re.findall(r'\d+', f))]
        
        print(f"Processing data in '{data_folder}'...")
        for f in files:
            path = os.path.join(data_folder, f)
            f_low = f.lower()
            
            sub = "others"
            if 'norma11' in f_low: sub = os.path.join("norma", "norma11")
            elif 'norma12' in f_low: sub = os.path.join("norma", "norma12")
            elif 'norma' in f_low: sub = "norma"
            elif any(v in f_low for v in ['u.txt', 'v.txt']): sub = "u_v"
            elif 'psi' in f_low: sub = "psi"
            elif 'tau' in f_low: sub = "tau"
            
            target = os.path.join(args.out, sub)
            
            if 'norma' in f_low:
                # Основной график жестких зон
                plot_field(path, f, target, special_mode=True, 
                           show_streamlines=True, tau_range=current_tau_range)
                # График с дугами
                plot_field(path, f, target, special_mode=True, 
                           show_streamlines=False, show_arcs=True, tau_range=current_tau_range)
            else:
                plot_field(path, f, target, special_mode=False, 
                           overlay_rigid=True, tau_range=current_tau_range)
        print("\nDone!")
