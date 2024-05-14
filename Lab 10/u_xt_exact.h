#ifndef U_XT_EXACT_H
#define U_XT_EXACT_H

#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<vector>
#include<complex>

using namespace std;

double u_xt_exact(double x, double t, double t0,  double freq, double sigma, double R, double G, double L, double C, 
			    double Rg, double Rl, double length, int number_nodes, int n_sum_terms);


#endif