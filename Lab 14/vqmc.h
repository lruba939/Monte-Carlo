#ifndef VQMC_H
#define VQMC_H

#include <vector>

double wav_f(double r, double a, double c);

double e_loc(double r, double a, double c);

double p_acc_fun(double r_i, double r_new, double a, double c);

std::vector<double> energy(double N, double dr, double r, double a, double c);

std::vector<double> histogram(double N, double dr, double r, double a, double c);



#endif