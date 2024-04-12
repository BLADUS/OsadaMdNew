#include "global_var.h"
#include "constants.h"
#include "params.h"
#include "start_cond.h"
#include <iostream>
#include <fstream>
#include <cmath> // Для использования функции sqrt() и pow()
#include <iomanip> // Для использования std::fixed и std::setprecision


//Выносим запись о информации кординат 
void writeCoordinateVector(std::ofstream& outputFile, int index) {
outputFile << "r" << index + 1 << " = (rx" << index + 1 << "; ry" 
<< index + 1 << "; rz" << index + 1 << ") = (" << coordx[index] << "; " 
<< coordy[index] << "; " << coordz[index] << ")" << std::endl;
}

//Выносим запись и вычисления о разности координат векторов 
void writeVectorDifference(std::ofstream& outputFile, int index1, int index2) {
    double diffx = coordx[index1] - coordx[index2];
    double diffy = coordy[index1] - coordy[index2];
    double diffz = coordz[index1] - coordz[index2];

    outputFile << "r" << index1 + 1 << index2 + 1 << " = (rx" << index1 + 1 << "; ry" << index1 + 1 << "; rz" << index1 + 1 << ") = (" << diffx << "; " << diffy << "; " << diffz << ")" << std::endl;
}

// Метод для вычисления и записи абсолютного значения векторов
void writeVectorAbsolute(std::ofstream& outputFile, int index1, int index2) {
    double diffx = coordx[index1] - coordx[index2];
    double diffy = coordy[index1] - coordy[index2];
    double diffz = coordz[index1] - coordz[index2];
    double distance = std::sqrt(diffx * diffx + diffy * diffy + diffz * diffz);

    outputFile << "r" << index1 + 1 << index2 + 1 << "_abs = " << distance << std::endl;
}

// Метод для записи информации о делении разности векторов на абсолютное значение
void writeVectorDivision(std::ofstream& outputFile, int index1, int index2) {
    double diffx = coordx[index1] - coordx[index2];
    double diffy = coordy[index1] - coordy[index2];
    double diffz = coordz[index1] - coordz[index2];
    double distance = std::sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
    
    // Запись разности векторов, деленной на абсолютное значение
    outputFile << "(rx" << index1 + 1 << index2 + 1 << "; ry" << index1 + 1 << index2 + 1 << "; rz" << index1 + 1 << index2 + 1 << ") / r" << index1 + 1 << index2 + 1 << "_abs = (" <<
        diffx / distance << "; "
        << diffy / distance << "; "
        << diffz / distance
        << ")" << std::endl;
}


void processAllVectorOperations(std::ofstream& outputFile, int numParticles) {
    for (int i = 0; i < numParticles; ++i) {
        for (int j = 0; j < numParticles; ++j) {
            if (i != j) {
                // Разность векторов
                writeVectorDifference(outputFile, i, j);
            }
        }
    }

    for (int i = 0; i < numParticles; ++i) {
        for (int j = 0; j < numParticles; ++j) {
            if (i != j) {
                // Абсолютное значение разности векторов
                writeVectorAbsolute(outputFile, i, j);
            }
        }
    }

    for (int i = 0; i < numParticles; ++i) {
        for (int j = 0; j < numParticles; ++j) {
            if (i != j) {
                // Деление векторов
                writeVectorDivision(outputFile, i, j);
            }
        }
    }
}



// Метод для записи информации о скорости частицы
void writeParticleVelocity(std::ofstream& outputFile, int index) {
    // Запись скорости частицы
    outputFile << "v" << index + 1 << " = (vx" << index + 1 << "; vy" << index + 1 << "; vz" << index + 1 << ") = (" << vx[index] << "; " << vy[index] << "; " << vz[index] << ")" << std::endl;
}

void calculate_forces(std::ofstream& outputFile) {
    // Проходим по парам частиц
    for (int i = 0; i < NUMBERPARTICLES - 1; ++i) {
        for (int j = i + 1; j < NUMBERPARTICLES; ++j) {
            // Вычисляем расстояние между частицами
            double distance = sqrt(pow(coordx[i] - coordx[j], 2) + pow(coordy[i] - coordy[j], 2) + pow(coordz[i] - coordz[j], 2));

            // Вычисляем силу между частицами
            double F = 24 * EPS / SIGMA * (2 * pow(SIGMA / distance, 13) - pow(SIGMA / distance, 7));

            // Рассчитываем составляющие силы по осям X, Y, Z
            double Fx_i = F * (coordx[i] - coordx[j]) / distance;
            double Fy_i = F * (coordy[i] - coordy[j]) / distance;
            double Fz_i = F * (coordz[i] - coordz[j]) / distance;

            double Fx_j = -Fx_i;
            double Fy_j = -Fy_i;
            double Fz_j = -Fz_i;

            // Записываем данные в файл
            outputFile << "F = " << F << "\n" // Запись силы между частицами
                       << "F" << i + 1 << " = (Fx" << i + 1 << "; Fy" << i + 1 << "; Fz" << i + 1 << ") = (" << Fx_i << "; " << Fy_i << "; " << Fz_i << ")\n" // Запись силы, действующей на первую частицу
                       << "F" << j + 1 << " = (Fx" << j + 1 << "; Fy" << j + 1 << "; Fz" << j + 1 << ") = (" << Fx_j << "; " << Fy_j << "; " << Fz_j << ")\n"; // Запись силы, действующей на вторую частицу
        }
    }
}

void calculate_potential_energy(std::ofstream& outputFile) {
    // Проходим по парам частиц
    for (int i = 0; i < NUMBERPARTICLES - 1; ++i) {
        for (int j = i + 1; j < NUMBERPARTICLES; ++j) {
            // Вычисляем расстояние между частицами
            double diffx = coordx[i] - coordx[j];
            double diffy = coordy[i] - coordy[j];
            double diffz = coordz[i] - coordz[j];
            double distance = std::sqrt(diffx * diffx + diffy * diffy + diffz * diffz);

            // Вычисляем потенциальную энергию между частицами
            double U12 = 4 * EPS * (std::pow(SIGMA / distance, 12) - std::pow(SIGMA / distance, 6));
            double U21 = U12; // Поскольку потенциальная энергия симметрична, U12 = U21

            // Записываем данные в файл
            outputFile << "U" << i + 1 << j + 1 << " = " << U12 << "\n"; // Запись потенциальной энергии между частицами
            outputFile << "U" << j + 1 << i + 1 << " = " << U21 << "\n"; // Запись потенциальной энергии между частицами
        }
    }
}

void write_to_file(){
    // Открытие файла для записи
    std::ofstream outputFile("Osada_MD_4.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл для записи\n";
        return;
    }
    
    // Установка точности для всего потока вывода
    outputFile << std::fixed << std::setprecision(8);
    // Запись в файл
    outputFile << "Step = 0" << std::endl;
    // Запись векторов координат в файл
    for (int i = 0; i < NUMBERPARTICLES; ++i) {
        writeCoordinateVector(outputFile, i);
    }

     processAllVectorOperations(outputFile, NUMBERPARTICLES);

    for (int i = 0; i < NUMBERPARTICLES; ++i) {
        writeParticleVelocity(outputFile, i);
    }
    
    calculate_potential_energy(outputFile);
    calculate_forces(outputFile);
    
    outputFile.close(); // Закрытие файла
}

void allocateMemory() {
    // Массив указателей на указатели
    double** arrays[] = {&coordx, &coordy, &coordz, &vx, &vy, &vz, &Fx, &Fy, &Fz};
    const int numArrays = sizeof(arrays) / sizeof(arrays[0]);

    // Выделение памяти для каждого массива
    for (int i = 0; i < numArrays; ++i) {
        *arrays[i] = (double*)malloc(NUMBERPARTICLES * sizeof(double));
        if (*arrays[i] == NULL) {
            std::cerr << "Ошибка: не удалось выделить память для массива " << i << std::endl;
            // Освобождение ранее выделенной памяти
            for (int j = 0; j < i; ++j) {
                free(*arrays[j]);
            }
            return;
        }
    }
}

// Функция для освобождения памяти
void freeMemory() {
    free(coordx);
    free(coordy);
    free(coordz);
    free(vx);
    free(vy);
    free(vz);
    free(Fx);
    free(Fy);
    free(Fz);
}

void simulation(){
    allocateMemory();
    start_cond_two_particles();
    write_to_file();
    freeMemory();
}

int main(){
    simulation();
}





