#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <math.h>
#include <vector>

double normal_rand();

double rand_in_range(double a, double b);

double g_function(double x, int i);

double fun_c_1(double x);

double fun_c_2(double x);

double fun_c_3(double x);

void basic_func(double a, double b, int k_start, int k_stop, int M, double (*function_C)(double), std::fstream *infile, std::fstream *log, std::vector<double>& con);

void opt_func(double a, double b, int k_start, int k_stop, int M, double **Nm_arr, int c, double (*function_C)(double), std::fstream *infile, std::fstream *log);

#endif