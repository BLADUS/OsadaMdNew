#include "global_var.h"
#include "constants.h"
#include "params.h"
#include "start_cond.h"
#include <iostream>
#include <fstream>
#include <cmath>   // Для использования функции sqrt() и pow()
#include <iomanip> // Для использования std::fixed и std::setprecision
#include <vector>

/////////////////////////////////////////////////UTIL/////////////////////////////////////////////////

// Выделение памяти
void allocateMemory()
{
    // Массив указателей на указатели
    double **arrays[] = {&coordx, &coordy, &coordz, &vx, &vy, &vz, &Fx, &Fy, &Fz};
    const int numArrays = sizeof(arrays) / sizeof(arrays[0]);

    // Выделение памяти для каждого массива
    for (int i = 0; i < numArrays; ++i)
    {
        *arrays[i] = (double *)malloc(NUMBERPARTICLES * sizeof(double));
        if (*arrays[i] == NULL)
        {
            std::cerr << "Ошибка: не удалось выделить память для массива " << i << std::endl;
            // Освобождение ранее выделенной памяти
            for (int j = 0; j < i; ++j)
            {
                free(*arrays[j]);
            }
            return;
        }
    }
}

// Функция для освобождения памяти
void freeMemory()
{
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
std::tuple<double, double, double> calculateCoordinateDifference(int index1, int index2)
{
    double diffx = coordx[index1] - coordx[index2];
    double diffy = coordy[index1] - coordy[index2];
    double diffz = coordz[index1] - coordz[index2];
    return std::make_tuple(diffx, diffy, diffz);
}

// Метод для вычисления и записи абсолютного значения векторов
double calculateVectorAbsolute(int index1, int index2)
{
    double diffx, diffy, diffz;
    std::tie(diffx, diffy, diffz) = calculateCoordinateDifference(index1, index2);
    double distance = std::sqrt(diffx * diffx + diffy * diffy + diffz * diffz);

    return distance;
}

// Метод для вычисления силы между двумя частицами
void calculateForce(int index1, int index2, double distance)
{
    if (distance <= RCUT)
    {
        F = (24 * EPS / distance) * (2 * std::pow((SIGMA / distance), 12) - std::pow((SIGMA / distance), 6));
    }
    else
    {
        F = 0;
    }
}

// Метод для вычисления вектора силы по осям X, Y, Z
void calculateForceVector(int index1, int index2, double distance)
{
    // Вычисляем силу между частицами
    calculateForce(index1, index2, distance);

    // Рассчитываем составляющие силы по осям X, Y, Z
    Fx[index1] = (F * (coordx[index1] - coordx[index2]) / distance) + 0.0;
    Fy[index1] = (F * (coordy[index1] - coordy[index2]) / distance) + 0.0;
    Fz[index1] = (F * (coordz[index1] - coordz[index2]) / distance) + 0.0;
    Fx[index2] = -Fx[index1];
    Fy[index2] = -Fy[index1];
    Fz[index2] = -Fz[index1];
}

// Метод для вычисления потенциальной энергии между двумя частицами
void calculatePotentialEnergy(int index1, int index2, double distance)
{
    // Вычисляем разность координат между частицами
    double diffx = coordx[index1] - coordx[index2];
    double diffy = coordy[index1] - coordy[index2];
    double diffz = coordz[index1] - coordz[index2];

    if (distance <= RCUT)
    {
        U = 4 * EPS * (std::pow(SIGMA / distance, 12) - std::pow(SIGMA / distance, 6)) - UCUT;
    }
    else
    {
        U = 0;
    }
}

void checkingBoundaryConditionCoordinates()
{
    for (int i = 0; i < NUMBERPARTICLES; i++)
    {
        // Проверка координаты x
        if (coordx[i] >= LX)
        {
            coordx[i] -= LX;
        }
        else if (coordx[i] < 0.0)
        {
            coordx[i] += LX;
        }

        // Проверка координаты y
        if (coordy[i] >= LY)
        {
            coordy[i] += -LY; // Отражение от верхней границы (можно изменить)
        }
        else if (coordy[i] < 0.0)
        {
            coordy[i] += LY; // Отражение от нижней границы (можно изменить)
        }

        // Проверка координаты z (аналогично x)
        if (coordz[i] >= LZ)
        {
            coordz[i] -= LZ;
        }
        else if (coordz[i] < 0.0)
        {
            coordz[i] += LZ;
        }
    }
}

void calculateVerleCoord()
{
    for (size_t i = 0; i < NUMBERPARTICLES; i++)
    {
        coordx[i] = coordx[i] + vx[i] * STEP + Fx[i] * std::pow(STEP, 2) / (2 * MASS);
        coordy[i] = coordy[i] + vy[i] * STEP + Fy[i] * std::pow(STEP, 2) / (2 * MASS);
        coordz[i] = coordz[i] + vz[i] * STEP + Fz[i] * std::pow(STEP, 2) / (2 * MASS);
    }
    checkingBoundaryConditionCoordinates();
}

void calculateVerleVelocity()
{
    for (size_t i = 0; i < NUMBERPARTICLES; i++)
    {
        vx[i] = vx[i] + Fx[i] / (2 * MASS) * STEP;
        vy[i] = vy[i] + Fy[i] / (2 * MASS) * STEP;
        vz[i] = vz[i] + Fz[i] / (2 * MASS) * STEP;
    }
}

/////////////////////////////////////////////////WRITE/////////////////////////////////////////////////
// Выносим запись о информации кординат
void writeCoordinateVector(std::ofstream &outputFile, int index)
{
    outputFile << "r" << index + 1 << " = (rx" << index + 1 << "; ry"
               << index + 1 << "; rz" << index + 1 << ") = (" << coordx[index] << "; "
               << coordy[index] << "; " << coordz[index] << ")" << std::endl;
}

// Запись разности координат векторов
void writeVectorDifference(std::ofstream &outputFile, int index1, int index2)
{
    double diffx, diffy, diffz;
    std::tie(diffx, diffy, diffz) = calculateCoordinateDifference(index1, index2);
    outputFile << "r" << index1 + 1 << index2 + 1 << " = (rx" << index1 + 1 << "; ry" << index1 + 1 << "; rz" << index1 + 1 << ") = (" << diffx << "; " << diffy << "; " << diffz << ")" << std::endl;
}

// Запись абсолютного значения векторов
void writeVectorAbsolute(std::ofstream &outputFile, int index1, int index2, double distance)
{
    outputFile << "r" << index1 + 1 << index2 + 1 << "_abs = " << distance << std::endl;
}

// Метод для записи информации о делении разности векторов на абсолютное значение
void writeVectorDivision(std::ofstream &outputFile, int index1, int index2, double distance)
{
    double diffx, diffy, diffz;
    std::tie(diffx, diffy, diffz) = calculateCoordinateDifference(index1, index2);
    // Запись разности векторов, деленной на абсолютное значение
    outputFile << "(rx" << index1 + 1 << index2 + 1 << "; ry" << index1 + 1 << index2 + 1 << "; rz" << index1 + 1 << index2 + 1 << ") / r" << index1 + 1 << index2 + 1 << "_abs = (" << diffx / distance << "; "
               << diffy / distance << "; "
               << diffz / distance
               << ")" << std::endl;
}

// Метод для записи информации о скорости частицы
void writeParticleVelocity(std::ofstream &outputFile, int index)
{
    // Запись скорости частицы
    outputFile << "v" << index + 1 << " = (vx" << index + 1 << "; vy" << index + 1 << "; vz" << index + 1 << ") = (" << vx[index] << "; " << vy[index] << "; " << vz[index] << ")" << std::endl;
}

// Метод для записи силы между частицами в файл
void writeForce(std::ofstream &outputFile, double distance)
{
    for (int i = 0; i < NUMBERPARTICLES - 1; ++i)
    {
        for (int j = i + 1; j < NUMBERPARTICLES; ++j)
        {
            calculateForce(i, j, distance);
            outputFile << "F" << i+1 << j+1 << " = " << F << "\n";
        }
    }
}

// Метод для записи вектора силы в файл
void writeForceVector(std::ofstream &outputFile)
{
    for (int i = 0; i < NUMBERPARTICLES - 1; ++i)
    {
        for (int j = i + 1; j < NUMBERPARTICLES; ++j)
        {
            outputFile << "F" << i + 1 << " = (Fx" << i + 1 << "; Fy" << i + 1
                       << "; Fz" << i + 1 << ") = (" << Fx[i] << "; " << Fy[i] << "; " << Fz[i] << ")\n";
        }
    }
}

// Метод для записи потенциальной энергии между частицами в файл
void writePotentialEnergy(std::ofstream &outputFile)
{
    // Проходим по парам частиц
    for (int i = 0; i < NUMBERPARTICLES - 1; ++i)
    {
        for (int j = i + 1; j < NUMBERPARTICLES; ++j)
        {
            // Записываем данные в файл
            outputFile << "U" << i + 1 << j + 1 << " = " << U << "\n"; // Запись потенциальной энергии между частицами
        }
    }
}

void write_data_step(std::ofstream &outputFile, int t, double distance)
{
    outputFile << "Step = " << t << std::endl;
    for (int i = 0; i < NUMBERPARTICLES; i++)
    {
        writeCoordinateVector(outputFile, i);
    }
    for (int i = 0; i < NUMBERPARTICLES - 1; ++i)
    {
        for (int j = i + 1; j < NUMBERPARTICLES; ++j)
        {
            writeVectorAbsolute(outputFile, i, j, distance);
        }
    }

    writePotentialEnergy(outputFile);
    writeForce(outputFile, distance);
    writeForceVector(outputFile);

    for (int i = 0; i < NUMBERPARTICLES; ++i)
    {
        writeParticleVelocity(outputFile, i);
    }

    outputFile << std::endl;
}

void calculate_step(std::ofstream &outputFile, int t)
{

    if (t != 0)
    {
        calculateVerleCoord();
        calculateVerleVelocity();
    }

    double distance = calculateVectorAbsolute(0, 1);
    calculateForce(0, 1, distance);
    calculatePotentialEnergy(0, 1, distance);

    // Рассчитываем составляющие силы по осям X, Y, Z
    Fx[0] = (F * (coordx[0] - coordx[1]) / distance) + 0.0;
    Fy[0] = (F * (coordy[0] - coordy[1]) / distance) + 0.0;
    Fz[0] = (F * (coordz[0] - coordz[1]) / distance) + 0.0;
    Fx[1] = -Fx[0];
    Fy[1] = -Fy[0];
    Fz[1] = -Fz[0];

    if (t != 0)
    {
        calculateVerleVelocity();
    }

    if (t <= 99 || t >= 900){
        write_data_step(outputFile, t, distance);
    }
}

void write_to_file()
{
    // Открытие файла для записи
    std::ofstream outputFile("Osada_MD_7.txt");
    if (!outputFile.is_open())
    {
        std::cerr << "Ошибка: не удалось открыть файл для записи\n";
        return;
    }

    // Установка точности для всего потока вывода
    outputFile << std::fixed << std::setprecision(8);
    for (int t = 0; t < NSTEPS; t++)
    {
        calculate_step(outputFile, t);
    }

    // Закрытие файла
    outputFile.close();
}

void simulation()
{
    allocateMemory();
    start_cond_two_particles();
    write_to_file();
    freeMemory();
}

int main()
{
    simulation();
}