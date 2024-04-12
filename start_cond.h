#pragma once
#include "global_var.h"

// Функция для задания ну двух частиц
void start_cond_two_particles() { 
	// Координаты
	coordx[0] = 14;	coordy[0] = 0.75;	coordz[0] = 0.5;
	coordx[1] = 14.5;	coordy[1] = 0.75;	coordz[1] = 0.5;
	// Скорости
	vx[0] = 1.0;	vy[0] = 1.0;	vz[0] = 0.0;
	vx[1] = -1.0;	vy[1] = 1.0;	vz[1] = 0.0;
}