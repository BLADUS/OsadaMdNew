#include "global_var.h"
#include "constants.h"
#include "params.h"
#include "start_cond.h"
#include <iostream>
#include <fstream>

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

void write_to_file(){
    // Открытие файла для записи
    std::ofstream outputFile("Osada_MD_2.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл для записи\n";
        return;
    }
    
    // Запись в файл
    outputFile.precision(8);
    outputFile << "\nNumber of Particles = " << NUMBERPARTICLES << std::endl;
    outputFile << "Number of Steps = " << NSTEPS << std::endl;
    outputFile << "Size of System on X = " << NUMCRIST_X << std::endl;
    outputFile << "Boltzmann Constant = " << K_B << std::endl;
    outputFile << "Volume = " << VOLUME << std::endl;
    outputFile.close(); // Закрытие файла
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
    write_to_file();
    freeMemory;
}

int main(){
    simulation();

    return 0;
}





