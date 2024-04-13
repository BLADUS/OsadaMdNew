#include "global_var.h"
#include "constants.h"
#include "params.h"
#include "start_cond.h"
#include <iostream>
#include <fstream>
#include <cmath> // Для использования функции sqrt() и pow()
#include <iomanip> // Для использования std::fixed и std::setprecision
#include <vector>

/////////////////////////////////////////////////UTIL/////////////////////////////////////////////////

// Выделение памяти
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

/////////////////////////////////////////////////CALCULATE/////////////////////////////////////////////////

// Метод для вычисления разности координат между двумя частицами
std::tuple<double, double, double> calculateCoordinateDifference(int index1, int index2) {
    double diffx = coordx[index1] - coordx[index2];
    double diffy = coordy[index1] - coordy[index2];
    double diffz = coordz[index1] - coordz[index2];
    return std::make_tuple(diffx, diffy, diffz);
}

// Метод для вычисления и записи абсолютного значения векторов
double calculateVectorAbsolute(int index1, int index2) {
    double diffx, diffy, diffz;
    std::tie(diffx, diffy, diffz) = calculateCoordinateDifference(index1, index2);
    double distance = std::sqrt(diffx * diffx + diffy * diffy + diffz * diffz);
    return distance;
}

/////////////////////////////////////////////////WRITE/////////////////////////////////////////////////
//Выносим запись о информации кординат 
void writeCoordinateVector(std::ofstream& outputFile, int index) {
outputFile << "r" << index + 1 << " = (rx" << index + 1 << "; ry" 
<< index + 1 << "; rz" << index + 1 << ") = (" << coordx[index] << "; " 
<< coordy[index] << "; " << coordz[index] << ")" << std::endl;
}

// Запись разности координат векторов 
void writeVectorDifference(std::ofstream& outputFile, int index1, int index2) {
    double diffx, diffy, diffz;
    std::tie(diffx, diffy, diffz) = calculateCoordinateDifference(index1, index2);

    outputFile << "r" << index1 + 1 << index2 + 1 << " = (rx" << index1 + 1 << "; ry" << index1 + 1 << "; rz" << index1 + 1 << ") = (" << diffx << "; " << diffy << "; " << diffz << ")" << std::endl;
}

// Запись абсолютного значения векторов
void writeVectorAbsolute(std::ofstream& outputFile, int index1, int index2) {
    double distance = calculateVectorAbsolute(index1, index2);
    outputFile << "r" << index1 + 1 << index2 + 1 << "_abs = " << distance << std::endl;
}


// Метод для записи информации о делении разности векторов на абсолютное значение
void writeVectorDivision(std::ofstream& outputFile, int index1, int index2) {
    double diffx, diffy, diffz;
    std::tie(diffx, diffy, diffz) = calculateCoordinateDifference(index1, index2);
    double distance = calculateVectorAbsolute(index1, index2);
    
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


void processAllVectorOperations(std::ofstream& outputFile, int numParticles) {
    // Проходим по всем парам частиц, но начинаем со второй частицы внешнего цикла, чтобы исключить дубликаты
    for (int i = 0; i < numParticles - 1; ++i) {
        for (int j = i + 1; j < numParticles; ++j) {
            // Записываем разность векторов и абсолютное значение разности
            writeVectorDifference(outputFile, i, j);
            writeVectorAbsolute(outputFile, i, j);
            writeVectorDivision(outputFile, i, j);
        }
    }
}

void write_to_file() {
    // Открытие файла для записи
    std::ofstream outputFile("Osada_MD_5.txt");
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

    // Вычисление и запись потенциальной энергии
    writePotentialEnergy(outputFile, NUMBERPARTICLES);
    
    // Вычисление и запись сил
    writeForces(outputFile);
    
    for (int i = 0; i < NUMBERPARTICLES; ++i) {
        writeParticleVelocity(outputFile, i);
    }
    
    // Закрытие файла
    outputFile.close();
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





