#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <bits/stdc++.h>
#include <iomanip>
#include <time.h>


double uni_rand();

void save_data(std::string outFileName, std::vector<std::vector<double>> data);

std::vector<std::vector<double>> load_data(std::string inFileName);

double R_0_f(double L, double C);

double zeta_f(double R_0, double R_g);

double Gamma_g_f(double R_0, double R_g);

double Gamma_l_f(double R_0, double R_l);

double c_f(double L, double C);

double mu_f(double G, double C);

double lambda_f(double R, double L, double G, double C);

double V_g(double t, double s);

double MC_line(int npaths, double x_start, double t_start, double c, double lambda, double mu,
             double zeta, double Gamma_g, double Gamma_l, double l);

#endif