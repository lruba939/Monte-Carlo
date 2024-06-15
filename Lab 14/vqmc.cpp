#include "utils.h"
#include "vqmc.h"

double wav_f(double r, double a, double c){
    double Psi = 0;
    Psi = (1+c*r)*exp(-a*r);

    return Psi;
}

double e_loc(double r, double a, double c){
    double e = 0;
    e = (-pow(a,2)*c*pow(r,2) + (-pow(a,2) + 4*a*c - 2*c)*r + 2*a - 2*c - 2) / (2*c*pow(r,2) + 2*r);

    return e;
}

double p_acc_fun(double r_i, double r_new, double a, double c){
    double p = 0;
    p = pow(r_new/r_i, 2) * pow(wav_f(r_new, a, c), 2) / pow(wav_f(r_i, a, c), 2);

    return p;
}

std::vector<double> energy(double N, double dr, double r, double a, double c){
    double r_i;
    r_i = r;
    std::vector<double> energy_stat(2);
    energy_stat[0] = 0; // mean
    energy_stat[1] = 0; // var

    for(int i=0; i<N; i++){
        double U1, U2;
        U1 = uni_rand();
        double r_new = r_i + dr*(2*U1 - 1);
        double p_acc = std::min(p_acc_fun(r_i, r_new, a, c), 1.);
        
        U2 = uni_rand();
        if(p_acc >= U2){
            r_i = r_new;
        }

        double energy = e_loc(r_i, a, c);
        energy_stat[0] += (double) energy;
        energy_stat[1] += (double) pow(energy,2);
    }

    energy_stat[0] = energy_stat[0]/N;
    energy_stat[1] = energy_stat[1]/N - pow(energy_stat[0]/N,2);

    return energy_stat;
}

std::vector<double> histogram(double N, double dr, double r, double a, double c){
    int M = 200;
    double r_max = 8.0;
    double delta_r = r_max/M;
    
    double r_i;
    r_i = r;
    std::vector<double> hist(M);

    for(int m = 0; m < M; m++){
        hist[m] = 0;
    }

    for(int i=0; i<N; i++){
        // std::cout<<i<<"\n";
        double U1, U2;
        U1 = uni_rand();
        double r_new = r_i + dr*(2*U1 - 1);
        double p_acc = std::min(p_acc_fun(r_i, r_new, a, c), 1.);
        
        U2 = uni_rand();
        if(p_acc >= U2){
            r_i = r_new;
        }

        if(r_i <= r_max){
            int k = (int) floor(r_i/delta_r);
            if(k < 200 && k >= 0){
                hist[k] += (double) 1/(N*delta_r);
            }
        }
    }

    return hist;
}