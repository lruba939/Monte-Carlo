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

std::vector<std::vector<double>> load_atoms(std::string inFileName);

double rho_q(double rho_max, double sigma_p, int i, int j, int n_x, int n_y, double Delta);

void relaksacja(const int n_x, const int n_y, double Delta, double epsilon, double omega, 
                double tol, double it_max, double rho_max, double sigma_p, double V_L,
                double V_T, double V_B);

void MC_V(const int n_x, const int n_y, double Delta, int N_chains_max, double n_length_max, double epsilon, 
                double rho_max, double sigma_p, double V_L, double V_T, double V_B);



#endif