import numpy as np
import matplotlib.pyplot as plt
import os
import re
import argparse
import gc

# ==========================================
# Настройки папок и визуализации
# ==========================================
DEFAULT_FOLDER = 'output_data'
DEFAULT_OUT_FOLDER = 'plots'
DEFAULT_TAU_S_RANGE = [0.2, 0.202] 
# ==========================================

def plot_convergence(filename, output_folder):
    """Отрисовка графика сходимости по итерациям из файла diff_max."""
    if not os.path.exists(filename):
        return
    
    try:
        # Загружаем данные. Файл может содержать один или несколько столбцов ошибок.
        data = np.loadtxt(filename)
        if data.size == 0: return

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')
        
        # Если в файле один столбец, приводим к 2D
        if len(data.shape) == 1:
            data = data.reshape(-1, 1)
            
        iterations = np.arange(1, data.shape[0] + 1)
        
        for i in range(data.shape[1]):
            label = f'Residual {i+1}' if data.shape[1] > 1 else 'Max Error'
            ax.semilogy(iterations, data[:, i], linewidth=2.5, label=label)
            
        ax.set_xlabel('Iteration', fontsize=18, fontweight='bold')
        ax.set_ylabel('Error Value (log)', fontsize=18, fontweight='bold')
        ax.tick_params(axis='both', which='major', labelsize=16)
        
        ax.grid(True, which="both", ls="--", alpha=0.5, color='gray')
        if data.shape[1] > 1:
            ax.legend(fontsize=14)
            
        plt.title('Algorithm Convergence', fontsize=22, pad=20, fontweight='bold')
        
        os.makedirs(output_folder, exist_ok=True)
        out_path = os.path.join(output_folder, 'convergence_log.png')
        plt.tight_layout()
        plt.savefig(out_path, dpi=300)
        plt.close()
        print(f"Convergence plot saved: {out_path}")
    except Exception as e:
        print(f"Error plotting convergence: {e}")

def load_data(filename):
    """Загрузка данных с оптимизацией памяти."""
    if not os.path.exists(filename):
        print(f"Файл не найден: {filename}")
        return None
    try:
        data = np.loadtxt(filename, dtype=np.float32)
        if len(data.shape) == 1:
            data = data.reshape(1, -1)
        return data
    except Exception:
        data = []
        with open(filename, 'r') as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        data.append([float(parts[0]), float(parts[1]), float(parts[2])])
                    except ValueError: continue
        return np.array(data, dtype=np.float32)

def get_best_colormap(z, special_mode=False):
    if special_mode:
        return plt.cm.plasma
    v_min, v_max = np.min(z), np.max(z)
    if v_min < 0 and v_max > 0: return 'RdBu_r'
    elif v_min >= 0: return 'magma'
    else: return 'plasma'

def plot_field(filename, title, output_folder='plots', special_mode=False, 
               show_contours=True, show_streamlines=True, overlay_rigid=False, tau_range=None,
               show_arcs=False):
    if tau_range is None:
        tau_range = DEFAULT_TAU_S_RANGE
        
    data = load_data(filename)
    if data is None or len(data) == 0:
        return

    x, y, z = data[:, 0], data[:, 1], data[:, 2]
    ux, uy = np.sort(np.unique(x)), np.sort(np.unique(y))
    nx, ny = len(ux), len(uy)

    # Проверка на достаточность данных для отрисовки (нужно минимум 2x2 или 3 точки для tricontour)
    if nx < 2 or ny < 2 or len(z) < 3:
        print(f"Пропуск {filename}: недостаточно данных для отрисовки ({nx}x{ny})")
        return

    # Для научного стиля используем квадратный формат и крупные шрифты (для всех графиков)
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

    cnt = None
    if len(z) != nx * ny:
        if special_mode:
            # Заливка жестких зон (Z < tau_s + eps) в научном стиле
            ax.tricontourf(x, y, z, levels=[0, tau_range[1]], colors=['#706b67'])
            ax.tricontour(x, y, z, levels=[tau_range[1]], colors='black', linewidths=1.5)
        else:
            limit = max(abs(z.min()), abs(z.max())) or 1.0
            vmin, vmax = (z.min(), z.max()) if cmap != 'RdBu_r' else (-limit, limit)
            cnt = ax.tricontourf(x, y, z, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        Z, X, Y = z.reshape((ny, nx)), *np.meshgrid(ux, uy)
        if special_mode:
            # Заливка жестких зон (Z < tau_s + eps)
            ax.contourf(X, Y, Z, levels=[0, tau_range[1]], colors=['#706b67'])
            # Тонкая черная граница жесткой зоны
            ax.contour(X, Y, Z, levels=[tau_range[1]], colors='black', linewidths=1.5)
        else:
            limit = max(abs(Z.min()), abs(Z.max())) or 1.0
            vmin, vmax = (Z.min(), Z.max()) if cmap != 'RdBu_r' else (-limit, limit)
            cnt = ax.contourf(X, Y, Z, levels=levels, cmap=cmap, vmin=vmin, vmax=vmax)
            if show_contours:
                ax.contour(X, Y, Z, levels=15, colors='black', linewidths=0.5, alpha=0.2)

    # Наложение линий тока (Серый цвет, отрисовка ПОВЕРХ жестких зон)
    if show_streamlines:
        dir_path = os.path.dirname(filename)
        match = re.search(r'(iter|iret)_?(\d+)', os.path.basename(filename))
        if match:
            prefix, num = match.group(1), match.group(2)
            psi_filename = None
            for pattern in [f"{prefix}{num}_psi.txt", f"{prefix}_{num}_psi.txt"]:
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

    # Наложение контуров жестких зон из файлов norma
    if overlay_rigid:
        dir_path = os.path.dirname(filename)
        match = re.search(r'(iter|iret)_?(\d+)', os.path.basename(filename))
        if match:
            prefix, num = match.group(1), match.group(2)
            # Ищем файлы типа iret101_norma11.txt, iret101_norma12.txt
            for n_suffix in ['norma11', 'norma12', 'norma']:
                for pattern in [f"{prefix}{num}_{n_suffix}.txt", f"{prefix}_{num}_{n_suffix}.txt"]:
                    n_path = os.path.join(dir_path, pattern)
                    if os.path.exists(n_path):
                        n_data = load_data(n_path)
                        if n_data is not None:
                            nx_n, ny_n = len(np.unique(n_data[:,0])), len(np.unique(n_data[:,1]))
                            if len(n_data) == nx_n * ny_n:
                                nz = n_data[:, 2].reshape((ny_n, nx_n))
                                nux, nuy = np.sort(np.unique(n_data[:,0])), np.sort(np.unique(n_data[:,1]))
                                ax.contour(nux, nuy, nz, levels=[tau_range[1]], 
                                           colors='black', linewidths=1.2, alpha=0.9)

    # Отрисовка теоретических дуг (если запрошено)
    if show_arcs:
        h_x = ux[1] - ux[0] if len(ux) > 1 else 0
        h_y = uy[1] - uy[0] if len(uy) > 1 else 0
        x0_bound, x1_bound = ux[0] - h_x/2, ux[-1] + h_x/2
        y0_bound, y1_bound = uy[0] - h_y/2, uy[-1] + h_y/2
        
        Lx = x1_bound - x0_bound
        # Радиус по новой формуле: Lx / (2 + pi)
        R = Lx / (2 + np.pi)
        
        t = np.linspace(0, np.pi/2, 100)
        # Отрисовка 4-х ВОГНУТЫХ дуг красным цветом (центры смещены внутрь)
        # Corner (x0, y0)
        ax.plot(x0_bound + R + R * np.cos(t + np.pi), y0_bound + R + R * np.sin(t + np.pi), 
                color='red', linewidth=4, linestyle='-', zorder=15)
        # Corner (x1, y0)
        ax.plot(x1_bound - R + R * np.cos(t + 3*np.pi/2), y0_bound + R + R * np.sin(t + 3*np.pi/2), 
                color='red', linewidth=4, linestyle='-', zorder=15)
        # Corner (x1, y1)
        ax.plot(x1_bound - R + R * np.cos(t), y1_bound - R + R * np.sin(t), 
                color='red', linewidth=4, linestyle='-', zorder=15)
        # Corner (x0, y1)
        ax.plot(x0_bound + R + R * np.cos(t + np.pi/2), y1_bound - R + R * np.sin(t + np.pi/2), 
                color='red', linewidth=4, linestyle='-', zorder=15)

    if cnt and not special_mode:
        fig.colorbar(cnt).ax.set_ylabel('Value', rotation=270, labelpad=15)
    
    ax.set_aspect('equal')
    plt.tight_layout()

    os.makedirs(output_folder, exist_ok=True)
    if show_arcs:
        suffix_base = "_rigid_arcs"
    else:
        suffix_base = "_yield_sci" if special_mode else "_full"
    out_name = os.path.join(output_folder, os.path.basename(filename).replace('.txt', f"{suffix_base}.png"))
    
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
    
    data_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), args.folder)
    if not os.path.exists(data_folder): data_folder = args.folder

    if os.path.exists(data_folder):
        # Отрисовка графика сходимости, если файл существует
        for conv_file in ['diff_max.txt', 'diff_max']:
            conv_path = os.path.join(data_folder, conv_file)
            if os.path.exists(conv_path):
                plot_convergence(conv_path, args.out)
                break

        # Исключаем файлы сходимости из общего цикла отрисовки полей
        files = sorted([f for f in os.listdir(data_folder) 
                        if f.endswith('.txt') and 'diff_max' not in f])
        if args.iter is not None:
            files = [f for f in files if any(int(n) == args.iter for n in re.findall(r'\d+', f))]
        
        print(f"Organizing plots into subfolders in '{args.out}'...")
        for f in files:
            path = os.path.join(data_folder, f)
            f_lower = f.lower()
            
            subfolder = ""
            is_norma = 'norma' in f_lower
            
            if 'norma11' in f_lower: subfolder = os.path.join("norma", "norma11")
            elif 'norma12' in f_lower: subfolder = os.path.join("norma", "norma12")
            elif 'norma' in f_lower: subfolder = "norma"
            elif any(v in f_lower for v in ['_u.txt', '_v.txt']):
                subfolder = "u_v"
            elif any(v in f_lower for v in ['_p.txt', 'r1']):
                subfolder = "p"
            elif '_psi.txt' in f_lower:
                subfolder = "psi"
            elif 'du_dx' in f_lower or 'dudx' in f_lower:
                subfolder = "du_dx"
            elif 'dv_dy' in f_lower or 'dvdy' in f_lower:
                subfolder = "dv_dy"
            elif 'ag_t' in f_lower:
                if 'ag_t11' in f_lower: subfolder = os.path.join("ag_t", "ag_t11")
                elif 'ag_t12' in f_lower: subfolder = os.path.join("ag_t", "ag_t12")
                else: subfolder = "ag_t"
            elif 'd12' in f_lower:
                subfolder = "D12"
                subfolder = "D12"
            elif 'tau' in f_lower:
                if 'tau11' in f_lower: subfolder = os.path.join("tau", "tau11")
                elif 'tau22' in f_lower: subfolder = os.path.join("tau", "tau22")
                elif 'tau12' in f_lower: subfolder = os.path.join("tau", "tau12")
                else: subfolder = "tau"
            elif 'gamma' in f_lower:
                if 'gamma11' in f_lower: subfolder = os.path.join("gamma", "gamma11")
                elif 'gamma22' in f_lower: subfolder = os.path.join("gamma", "gamma22")
                elif 'gamma12' in f_lower: subfolder = os.path.join("gamma", "gamma12")
                else: subfolder = "gamma"
            elif 'du_dy' in f_lower:
                subfolder = "du_dy"
            elif 'dv_dx' in f_lower:
                subfolder = "dv_dx"
            else:
                subfolder = "others"
                
            target_dir = os.path.join(args.out, subfolder)
            
            if is_norma:
                # 1. Качественный график с линиями тока (psi) и жесткими зонами
                plot_field(path, f, target_dir, special_mode=True, 
                           show_streamlines=True, tau_range=current_tau_range)
                # 2. ОТДЕЛЬНАЯ КАРТИНКА: Жесткие зоны + теоретические дуги (красные)
                plot_field(path, f, target_dir, special_mode=True, 
                           show_streamlines=False, show_arcs=True, tau_range=current_tau_range)
            else:
                # Для всех остальных полей: отключаем psi, но накладываем границы жестких зон
                plot_field(path, f, target_dir, special_mode=False, 
                           show_streamlines=False, overlay_rigid=True, tau_range=current_tau_range)
        print("\nDone! All plots are organized.")
