#include "functions.h"



int main(){

    double dt = 0.1, tmax = 100;
    double D=1;
    double Exy = 0;
      
    std::fstream wiener_calc;
    wiener_calc.open("wiener.dat", std::fstream::out);

    for(int N=pow(10,2); N<=pow(10,5); N*=10){
        wiener(N, tmax, dt, Exy, D, &wiener_calc);
    }
    wiener_calc.close();
    

    return 0;

}