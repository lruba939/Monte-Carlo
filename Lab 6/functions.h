#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "particle_translation.h"
#include <fstream>
#include <iostream>


double uni_rand();

void normal_rand(double& X, double& Y, double std);

void wiener(int N, double tmax, double dt, double D, std::fstream *pos, std::fstream *stat);

void dyf_and_abs(int N, double tmax, double dt, double xs, double ys, double D, double omega,
                 double xr, double yr, double Rr, double xa, double ya, double Ra,
                 std::fstream *pos, std::fstream *stat);

#endif