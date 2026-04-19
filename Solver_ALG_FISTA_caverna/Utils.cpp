#include "Utils.h"

Timer::Timer() : start_time(std::chrono::high_resolution_clock::now()) {}

void Timer::reset() {
    start_time = std::chrono::high_resolution_clock::now();
}

double Timer::elapsed() const {
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    return diff.count();
}

void reading(const std::string& filename, std::vector<std::vector<double>>& data1)
{
    // 1. Открываем файл для чтения
    std::ifstream file(filename);

    // Проверка, удалось ли открыть файл
    if (!file.is_open()) {
        std::cerr << "Error opening file for reading: " << filename << std::endl;
        return;
    }

    int i, j;
    double value;

    // 2. Читаем данные пока файл не закончится (берем i, j и value из каждой строки)
    while (file >> i >> j >> value)
    {
        // Проверка границ массива (чтобы не выйти за пределы вектора)
        if (j >= 0 && j < (int)data1.size()) {
            if (i >= 0 && i < (int)data1[j].size()) {
                // Заполняем массив: в вашем случае индекс j — строки, i — столбцы
                data1[j][i] = value;
            }
        }
    }

    // 3. Закрываем файл
    file.close();
}

void save_iteration(const std::string& filename, const std::vector<std::vector<double>>& data1, int nx, int ny) 
{

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for saving: " << filename << std::endl;
        return;
    }
    file << std::fixed << std::setprecision(8);

    
        for (int j = 1; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i) 
                file << i << " " << j << " " << data1[j][i] << std::endl;

            //file << std::endl;
        } 
    file.close();
}

void save_field(const std::string& filename, const std::vector<std::vector<double>>& data1, int nx, int ny, double hx, double hy, int boll) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file for saving: " << filename << std::endl;
        return;
    }
    file << std::fixed << std::setprecision(8);

    int Nx1 = nx + 2;
    int Ny1 = ny + 2;
    double hx_2 = hx / 2.0;
    double hy_2 = hy / 2.0;
    double Lx1 = nx * hx;
    double Ly1 = ny * hy;
    double x, y;

    if (boll == 1) { // u (vertical faces)
        for (int j = 1; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                x = (i)*hx;
                y = (j)*hy - hy_2;
                file << x << " " << y << " " << data1[j][i] << std::endl;
            }
            file << std::endl;
        }
    }
    else if (boll == 2) { // v (horizontal faces)
        for (int j = 0; j <= ny; ++j) {
            for (int i = 1; i <= nx; ++i) {
                x = (i)*hx - hx_2;
                y = (j)*hy;
                file << x << " " << y << " " << data1[j][i] << std::endl;
            }
            file << std::endl;
        }
    }
    else if (boll == 3) { // centers (p, tau, etc.)
        for (int j = 1; j <= ny; ++j) {
            for (int i = 1; i <= nx; ++i) {
                x = (i)*hx - hx_2;
                y = (j)*hy - hy_2;
                file << x << " " << y << " " << data1[j][i] << std::endl;
            }
            file << std::endl;
        }
    }
    else if (boll == 4) { // nodes (psi)
        for (int j = 0; j <= ny; ++j) {
            for (int i = 0; i <= nx; ++i) {
                x = (i)*hx;
                y = (j)*hy;
                file << x << " " << y << " " << data1[j][i] << std::endl;
            }
            file << std::endl;
        }
    }
    else if (boll == 5) {

        file << 0.0 << " " << 0.0 << " " << data1[0][0] << std::endl;
        for (int i = 1; i <= nx; ++i)
        {
            x = i * hx - hx / 2.0; y = 0;
            file << x << " " << y << " " << data1[0][i] << std::endl;
        }
        file << Lx1 << " " << 0.0 << " " << data1[0][nx + 1] << std::endl;
        for (int j = 1; j <= ny; ++j) {
            x = 0; y = j * hy - hy / 2.0;
            file << x << " " << y << " " << data1[j][0] << std::endl;

            for (int i = 1; i <= nx; ++i) {
                x = i * hx - hx / 2.0; y = j * hy - hy / 2.0;
                file << x << " " << y << " " << data1[j][i] << std::endl;
            }
            x = Lx1; y = j * hy - hy / 2.0;
            file << x << " " << y << " " << data1[j][nx + 1] << std::endl;
            file << std::endl;
        }
        file << 0.0 << " " << Ly1 << " " << data1[ny + 1][0] << std::endl;
        for (int i = 1; i <= nx; ++i) {
            
            x = i * hx - hx / 2.0; y = Ly1;
            file << x << " " << y << " " << data1[ny+1][i] << std::endl;
        }
        // Corners
        file << Lx1 << " " << Ly1 << " " << data1[ny + 1][nx + 1] << std::endl;
    }

    file.close();
}
