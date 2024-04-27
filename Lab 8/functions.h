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


double uni_rand();

void save_data(std::string outFileName, std::vector<std::vector<double>> data);

std::vector<std::vector<double>> load_atoms(std::string inFileName);

double cut_function(double r, double R1, double R2);

double B_val(std::vector<std::vector<double>> atoms, int i, int j, double r_ij[3], double r_ij_len, double R1, double R2, double a0, double c0, double d0, double delta);

double Brenner_tot_energy(std::vector<std::vector<double>> atoms_xyz);

double Brenner_i_energy(std::vector<std::vector<double>> atoms_xyz, int i);

std::vector<double> pair_correlation_function(std::vector<std::vector<double>> atoms, int M);

std::vector<std::vector<double>> rand_atoms_pos(int n, double r);

std::vector<std::vector<double>> atom_pos_change(std::vector<std::vector<double>> atoms_rad, int i, double beta, double w_r, double w_phi, double w_theta);

std::vector<std::vector<double>> atoms_global_pos_change(std::vector<std::vector<double>> atoms_rad, double beta, double W_all, double& V_tot);

void SA(int n, int it_max, int M, double w_r, double w_phi, double w_theta, double W_all, double r_start, double beta_min, double beta_max, double p, std::ofstream *energyInFile, std::ofstream *posInFile, std::ofstream *pcfInFile);

std::vector<std::vector<double>> rad_to_xyz(std::vector<std::vector<double>> atoms_rad);

double beta_parameter(int it, int it_max, double beta_min, double beta_max, double p);

#endif