#ifndef FUNCTIONS_2_H
#define FUNCTIONS_2_H

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

double cut_function(double r, double R1, double R2);

double B_val(std::vector<std::vector<double>>& atoms, int i, int j, double r_ij[3], double r_ij_len, double R1, double R2, double a0, double c0, double d0, double delta, bool correction);

double Brenner_i_energy(std::vector<std::vector<double>>& atoms, int i, bool correction);

std::vector<double> pair_correlation_function(std::vector<std::vector<double>>& atoms, int M, double& r_mean);

void rand_atoms_pos(std::vector<std::vector<double>>& atoms, const int n, double r);

void atom_pos_change(std::vector<std::vector<double>>& atoms, int i, double beta, double w_r, double w_phi, double w_theta, bool correction);

void atoms_global_pos_change(std::vector<std::vector<double>>& atoms , const int n, double beta, double W_all, double& V_tot, bool correction);

void rad_to_xyz(std::vector<std::vector<double>>& atoms, int i);

double beta_parameter(int it, int it_max, double beta_min, double beta_max, double p);

double SA(const int n, const int it_max, int M, double w_r, double w_phi, double w_theta, double W_all, double r_start, double beta_min, double beta_max, double p,
        bool correction, std::ofstream *energyInFile, std::ofstream *posInFile, std::ofstream *pcfInFile, std::ofstream *avogadroInFile);


#endif