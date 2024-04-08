#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "particle_translation.h"
#include <fstream>
#include <iostream>


double uni_rand();

void normal_rand(double& X, double& Y, double mean, double std);

void wiener(int N, double tmax, double dt, double Exy, double D, std::fstream *pos, std::fstream *stat);

#endif