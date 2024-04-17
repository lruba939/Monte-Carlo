#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <fstream>
#include <iostream>
#include <cmath>


double uni_rand();

void normal_rand(double& X, double& Y, double mean, double std);

void Gillespie(double k1, double k2, double k3, double k4, double x1start, double x2start, double x3start, int tmax, int N, int Pmax, std::fstream *stat, std::fstream *hist);


#endif