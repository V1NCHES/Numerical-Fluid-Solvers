import numpy as np
import matplotlib.pyplot as plt
import os
import re
import argparse
import gc

# ==========================================
# Настройки по умолчанию
# ==========================================
DEFAULT_FOLDER = 'output_data_boost'
DEFAULT_OUT_FOLDER = 'plots'
DEFAULT_TAU_S_RANGE = [0.1*np.sqrt(2), 0.1*np.sqrt(2)+0.01] 
# ==========================================

# L-shape geometry (must match C++ config)
NX1, NX2, NY1, NY2 = 100, 200, 100, 200
H = 0.01


def load_data(filename):
    """Загрузить файл в массив numpy. Пустые строки используются как разделители строк."""
    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        return None
    try:
        data = np.loadtxt(filename, dtype=np.float64)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return None


def is_in_L(x, y, nx1=NX1, nx2=NX2, ny1=NY1, ny2=NY2, h=H, eps=1e-9):
    """Проверка, попадает ли точка (x,y) в L-образную область."""
    Lx1, Lx2 = nx1 * h, nx2 * h
    Ly1, Ly2 = ny1 * h, ny2 * h
    if 0 - eps <= x <= Lx2 + eps and 0 - eps <= y <= Ly1 + eps:
        return True
    if 0 - eps <= x <= Lx1 + eps and Ly1 - eps <= y <= Ly2 + eps:
        return True
    return False


def plot_field(filename, out_folder, title=None, cmap='viridis', levels=50, special_mode=False, tau_range=None):
    """Загрузить файл и нарисовать поле на L-образной области."""
    if tau_range is None: tau_range = DEFAULT_TAU_S_RANGE
    
    data = load_data(filename)
    if data is None or len(data) == 0:
        return

    x_raw, y_raw, z_raw = data[:, 0], data[:, 1], data[:, 2]

    # Определяем уникальные координаты сетки
    ux = np.sort(np.unique(x_raw))
    uy = np.sort(np.unique(y_raw))
    
    # Создаем 2D сетку
    X, Y = np.meshgrid(ux, uy)
    Z = np.full(X.shape, np.nan) # Инициализируем пустотой

    # Заполняем сетку данными
    # Создаем словарь для быстрого поиска: (x, y) -> z
    lookup = {(xi, yi): zi for xi, yi, zi in zip(x_raw, y_raw, z_raw)}
    
    for i in range(len(uy)):
        for j in range(len(ux)):
            xi, yi = ux[j], uy[i]
            if is_in_L(xi, yi):
                Z[i, j] = lookup.get((xi, yi), np.nan)

    if np.all(np.isnan(Z)):
        print(f"No valid points in L-domain for {filename}")
        return

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('white')

    if special_mode:
        # Заливка жестких зон (Z < tau_s) 
        # Используем pcolormesh для дискретного отображения зон
        Z_mask = np.where(Z < tau_range[1], 1.0, 0.0)
        Z_mask[np.isnan(Z)] = np.nan
        
        cnt = ax.pcolormesh(X, Y, Z_mask, cmap='Greys', shading='auto', vmin=0, vmax=2)
        # Добавляем контур границы жесткой зоны
        ax.contour(X, Y, Z, levels=[tau_range[1]], colors='red', linewidths=1.5)
    else:
        # Используем pcolormesh для качественного отображения поля
        cnt = ax.pcolormesh(X, Y, Z, cmap=cmap, shading='auto')
        ax.contour(X, Y, Z, levels=15, colors='black', linewidths=0.4, alpha=0.3)
        cb = fig.colorbar(cnt, ax=ax, fraction=0.03, pad=0.04)
        cb.ax.set_ylabel('Value', rotation=270, labelpad=15, fontsize=12)

    # Draw L-shape boundary
    Lx1, Lx2 = NX1 * H, NX2 * H
    Ly1, Ly2 = NY1 * H, NY2 * H
    boundary_coords = [(0,0),(Lx2,0),(Lx2,Ly1),(Lx1,Ly1),(Lx1,Ly2),(0,Ly2),(0,0)]
    bx, by = zip(*boundary_coords)
    ax.plot(bx, by, color='black', linewidth=2.5)

    ftitle = title or os.path.basename(filename).replace('.txt', '')
    ax.set_title(ftitle, fontsize=16, fontweight='bold', pad=12)
    ax.set_xlabel('x', fontsize=13)
    ax.set_ylabel('y', fontsize=13)
    ax.set_aspect('equal')
    
    # Устанавливаем лимиты осей чуть шире области
    ax.set_xlim(-H, Lx2 + H)
    ax.set_ylim(-H, Ly2 + H)
    
    plt.tight_layout()

    os.makedirs(out_folder, exist_ok=True)
    out_name = os.path.join(out_folder, os.path.basename(filename).replace('.txt', '_clean.png'))
    plt.savefig(out_name, dpi=200)
    plt.close(fig)
    print(f"Saved: {out_name}")
    gc.collect()


def plot_convergence(filename, out_folder):
    """Plot convergence from diff_max.txt / similar files."""
    if not os.path.exists(filename):
        return
    try:
        data = np.loadtxt(filename)
        if data.size == 0:
            return
        if data.ndim == 1:
            data = data.reshape(-1, 1)
        fig, ax = plt.subplots(figsize=(10, 6))
        iters = np.arange(1, data.shape[0] + 1)
        for col in range(data.shape[1]):
            label = f'Residual {col+1}' if data.shape[1] > 1 else 'Max Error'
            ax.semilogy(iters, data[:, col], linewidth=2, label=label)
        ax.set_xlabel('Iteration', fontsize=14)
        ax.set_ylabel('Error (log)', fontsize=14)
        ax.grid(True, which='both', ls='--', alpha=0.5)
        if data.shape[1] > 1:
            ax.legend(fontsize=12)
        ax.set_title('Convergence', fontsize=16, fontweight='bold')
        plt.tight_layout()
        os.makedirs(out_folder, exist_ok=True)
        out_path = os.path.join(out_folder, 'convergence.png')
        plt.savefig(out_path, dpi=200)
        plt.close()
        print(f"Convergence plot: {out_path}")
    except Exception as e:
        print(f"Error plotting convergence: {e}")


# ==========================================
# Field type detection by filename
# ==========================================
FIELD_CMAPS = {
    'u':      ('viridis',   'U Velocity'),
    'v':      ('plasma',    'V Velocity'),
    'p':      ('coolwarm',  'Pressure'),
    'psi':    ('RdBu_r',    'Stream Function'),
    'tau':    ('hot',       'Tau'),
    'gamma':  ('YlOrRd',    'Gamma'),
    'norma':  ('jet',       'Normal Stress'),
    'ag_t':   ('inferno',   'Augmented Tensor'),
}

def detect_cmap(fname):
    fl = fname.lower()
    for key, (cmap, title) in FIELD_CMAPS.items():
        if key in fl:
            return cmap, title
    return 'viridis', 'Field'


def get_subfolder(fname):
    fl = fname.lower()
    for key in ['norma11','norma12','norma','ag_t11','ag_t12','ag_t',
                'tau11','tau22','tau12','tau',
                'gamma11','gamma22','gamma12','gamma',
                'psi','_u.txt','_v.txt','_p.txt']:
        if key in fl:
            clean = key.replace('.txt','').strip('_')
            return clean
    return 'others'


# ==========================================
# Main
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot L-shape solver output fields')
    parser.add_argument('--folder', type=str, default=DEFAULT_FOLDER)
    parser.add_argument('--out', type=str, default=DEFAULT_OUT_FOLDER)
    parser.add_argument('--tau_min', type=float, default=DEFAULT_TAU_S_RANGE[0])
    parser.add_argument('--tau_max', type=float, default=DEFAULT_TAU_S_RANGE[1])
    args = parser.parse_args()
    
    current_tau_range = [args.tau_min, args.tau_max, DEFAULT_TAU_S_RANGE[1]]

    data_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), args.folder)
    if not os.path.exists(data_folder):
        data_folder = args.folder

    if not os.path.exists(data_folder):
        print(f"Folder not found: {data_folder}")
        exit(1)

    # Convergence
    for conv_file in ['diff_max.txt', 'diff_max']:
        conv_path = os.path.join(data_folder, conv_file)
        if os.path.exists(conv_path):
            plot_convergence(conv_path, args.out)
            break

    # All .txt fields
    files = sorted([f for f in os.listdir(data_folder)
                    if f.endswith('.txt') and 'diff_max' not in f])
    print(f"Found {len(files)} field files in '{data_folder}'")

    for f in files:
        path = os.path.join(data_folder, f)
        cmap, title = detect_cmap(f)
        subfolder = os.path.join(args.out, get_subfolder(f))
        
        is_norma = 'norma' in f.lower()
        plot_field(path, subfolder, title=f"{title}: {f}", cmap=cmap, special_mode=is_norma, tau_range=current_tau_range)

    print("\nDone! All plots saved to:", args.out)
