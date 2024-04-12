#pragma once
// Координаты ,массивы от числа частиц, 
double* coordx;  // *  для выделения памяти (malloc)
double* coordy;
double* coordz;
// Скорости (массивы от числа частиц, i - номер частицы) 
double* vx;
double* vy;
double* vz;
// Потенциальная, кин., кин. тепл., полная, внутр. энергии(скаляры)
double U;    // потенциальная энергия
double Ekin; // кин. энергия
double Eterm;// кин. тепл. энергия
double E;    // полная энергия
double Eint; // внутр. энергия
// Силы взаимодействия (массивы от числа частиц и вспомогательный скаляр)
double* Fx;
double* Fy;
double* Fz;
double F;
// Энергии (ср. на 1 частицу)
double U1;
double Ekin1;
double Eterm1;
double E1;
double Eint1;
// Энергии (ср. на 1 частицу и на шаг)
double U1_av;
double Ekin1_av;
double Eterm1_av;
double E1_av;
double Eint1_av;
// Темапература
double T;
double T_av;
// Давление 
double P;
double P1;
double P2;
double P_av;
// Импульс
double p;
double px;
double py;
double pz;