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
    outputFile << "r" << index + 1 << " = (" << coordx[index] << "; " << coordy[index] << "; " << coordz[index] << ")" << std::endl;
}

//Выносим запись и вычисления о разности координат векторов 
void writeVectorDifference(std::ofstream& outputFile, int index1, int index2) {
    double diffx = coordx[index1] - coordx[index2];
    double diffy = coordy[index1] - coordy[index2];
    double diffz = coordz[index1] - coordz[index2];

    outputFile << "r" << index1 + 1 << index2 + 1 << " = (" << diffx << "; " << diffy << "; " << diffz << ")" << std::endl;
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

// Метод для записи информации о скорости частицы
void writeParticleVelocity(std::ofstream& outputFile, int index) {
    // Запись скорости частицы
    outputFile << "v" << index + 1 << " = (vx" << index + 1 << "; vy" << index + 1 << "; vz" << index + 1 << ") = (" << vx[index] << "; " << vy[index] << "; " << vz[index] << ")" << std::endl;
}

void write_to_file(){
    // Открытие файла для записи
    std::ofstream outputFile("Osada_MD_3.txt");
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

    writeVectorDifference(outputFile, 0, 1);
    writeVectorDifference(outputFile, 1, 0);

    writeVectorAbsolute(outputFile, 0, 1);
    writeVectorAbsolute(outputFile, 0, 1);

    writeVectorDivision(outputFile, 0, 1);
    writeVectorDivision(outputFile, 1, 0);

    for (int i = 0; i < NUMBERPARTICLES; ++i) {
        writeParticleVelocity(outputFile, i);
    }
    
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





