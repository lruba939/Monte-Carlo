#include "functions.h"



int main(){

    double dt = 0.1, tmax = 100;
    double D=1;
    double Exy = 0;
      
    std::fstream wiener_calc, wiener_stat;
    wiener_calc.open("wiener_pos.dat", std::fstream::out);
    wiener_stat.open("wiener_stat.dat", std::fstream::out);
        wiener_stat<<"N t E_X E_Y E_XY D_XX D_YY D_XY\n";

    for(int N=pow(10,2); N<=pow(10,5); N*=10){
        wiener(N, tmax, dt, Exy, D, &wiener_calc, &wiener_stat);
    }
    wiener_calc.close();
    

    return 0;

}