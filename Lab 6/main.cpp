#include "functions.h"

#include <cstdlib>
#include <cstdio>
#include <ctime>

int main(){

    double dt = 0.1, tmax = 100;
    double D=1;
    double xs = 0, ys = 0;
      
    std::fstream wiener_calc, wiener_stat;
    wiener_calc.open("wiener_pos.dat", std::fstream::out);
    wiener_stat.open("wiener_stat.dat", std::fstream::out);
        wiener_stat<<"N t E_X E_Y E_XY D_XX D_YY D_XY\n";

    for(int N=pow(10,2); N<=pow(10,5); N*=10){
        wiener(N, tmax, dt, D, &wiener_calc, &wiener_stat);
    }
    wiener_calc.close();
    wiener_stat.close();
    
    // D = 1;
    // double xr = 0, yr = 0, Rr = 5;
    // double xa = 3, ya = 0, Ra[2] = {0.1, 0.5};
    // xs = -4.5; ys = 0;
    // int Nmax = pow(10,4);
    // tmax = 1000; dt = 0.1;
    // double omega[3] = {10, 50, 100};

    // std::fstream da_pos, da_stat;
    // da_pos.open("omega2_Ra_1_pos.dat", std::fstream::out);
    // da_stat.open("omega2_Ra_1_stat.dat", std::fstream::out);
    //     dyf_and_abs(Nmax, tmax, dt, xs, ys, D, omega[2], xr, yr, Rr, xa, ya, Ra[1], &da_pos, &da_stat);
    // da_pos.close();
    // da_stat.close();

    return 0;

}