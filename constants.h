#pragma once
#include <cmath>
#include "params.h"
const double MASS= 66.335; // масса одной частицы (атом аргона)
const double K_B = 1.380648528; // постоянная Больцмана

// параметры потенциала Л.-Дж. для аргона
const double EPS = 1.712; // параметр энергии (глубина потенциальной ямы)
const double SIGMA = 0.3418; // параметр длины взаимодействия
const double RCUT = 2.5 * SIGMA; // радиус обрезания потенциала
const double RCUT2 = RCUT * RCUT;
const double UCUT = 4 * EPS * (pow((SIGMA / RCUT), 12) - pow((SIGMA / RCUT), 6));

const double ACRIST = 1.0; // длина ребра эл. ячейки (зависит от термодинамического состояния и модели)
// размеры системы по осям координат
const double LX = NUMCRIST_X * ACRIST;
const double LY = NUMCRIST_Y * ACRIST;
const double LZ = NUMCRIST_Z * ACRIST;
const double VOLUME = LX * LY * LZ; // объем системы