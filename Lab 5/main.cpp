#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <math.h>


double normal_rand(){
    return (rand()/(double) RAND_MAX);
}

double rand_in_range(double a, double b){
    return a + (b - a)*normal_rand();
}

double g_function(double x, int i){
    switch (i){
        case 1:
            return 1. + tanh(x);

        case 2:
            return 1./(1+pow(x,2));
        
        case 3:
            return pow(cos(M_PI*x),10);
        
        default:
            std::cout<<"No matching function.";
            return 0;
    }
}

double ranges[3][3] = {
    {-3, 3}, // całka C1
    {0, 10}, // całka C2
    {0, 1} // całka C3
};


int main(){

    double xi = 0; double result = 0;

    // Metoda podstawowa

    std::fstream podst;
    podst.open("podst.dat", std::fstream::out);

    // Pętla odpowiedzialna za wybór funkcji
    for(int c=1; c<=3; c++){
        
        double a = ranges[c-1][0];
        double b = ranges[c-1][1];

        //Pętla odpowiedzialna wybór ilości losowań
        for(int k=2; k <= 5; k++){
            
            int N = pow(10,k);
            double stat_values[5] = {0,0,0,0,0}; // {mean, second moment, var, std, error}

            // Pętla odpowiedzialna za losowanie
            for(int i=0; i<N; i++){
                xi = rand_in_range(a, b);
                result = (b - a)*g_function(xi, c);
                stat_values[0] += result/N;
                stat_values[1] += pow(result, 2)/N;
            }
            stat_values[2] = (stat_values[1] - pow(stat_values[0], 2))/N;
            stat_values[3] = sqrt(stat_values[2]);
            stat_values[4] = stat_values[3]/stat_values[0]*100;

            podst<<c<<" "<<k<<" "<<stat_values[0]<<" "<<stat_values[1]<<" "<<stat_values[2]<<" "<<stat_values[3]<<" "<<stat_values[4]<<"\n";
        }
    }
    podst.close();


    // Metoda losowania systematycznego (warstwowe nieoptymalne)
    std::fstream sys;
    sys.open("sys.dat", std::fstream::out);

    int M = 10;
    double pm = 1./M;

    for(int c=1; c<=3; c++){
        double a = ranges[c-1][0];
        double b = ranges[c-1][1];
        double delta_x = (b - a)/M;

        for(int k=2; k <= 5; k++){
            int N = pow(10,k);
            double stat_values[5] = {0,0,0,0,0}; // {mean, second moment, var, std, error}
            double Nm = N*pm;

            for(int m=0; m<M; m++){
                double x_start = a + delta_x*m;
                double x_stop = x_start + delta_x;
                // std::cout<<"x_start = "<<x_start<<"\n";

                for(int i=0; i<Nm; i++){
                    xi = rand_in_range(x_start, x_stop);
                    result = (b - a)*g_function(xi, c);
                    stat_values[0] += result/Nm;
                    stat_values[1] += pow(result, 2)/Nm;
                    stat_values[2] += (stat_values[1] - pow(stat_values[0], 2));
                }
            }
            stat_values[0] = stat_values[0]*pm;
            stat_values[1] = stat_values[1]*pm;
            stat_values[2] = (stat_values[1] - pow(stat_values[0], 2))/Nm*pm*pm;
            stat_values[3] = sqrt(stat_values[2]);
            stat_values[4] = stat_values[3]/stat_values[0]*100;

            sys<<c<<" "<<k<<" "<<stat_values[0]<<" "<<stat_values[1]<<" "<<stat_values[2]<<" "<<stat_values[3]<<" "<<stat_values[4]<<"\n";
        }
    }
    sys.close();


    return 0;
}