#pragma once
const int NSTEPS = 1000; // число шагов
const int LASTSTEP = NSTEPS - 1; // последний шаг

// число элементарных ячеек (кристаллов) по осям координат
const int NUMCRIST_X = 2; 
const int NUMCRIST_Y = 2; 
const int NUMCRIST_Z = 2;

// число частиц для примитивной решетки с учетом ПГУ (периодические граничные условия)
//NUMBERPARTICLES = NUMCRIST_X * NUMCRIST_Y * NUMCRIST_Z;
//NUMBERPARTICLES = <значение>;

const int NUMBERPARTICLES = 2;
const double STEP = 0.002; // шаг интегрирования разностной схемы