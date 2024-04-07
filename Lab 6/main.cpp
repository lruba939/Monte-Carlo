#include "particle_translation.h"
#include <random>
#include <fstream>


int main(){

    int N = 10000;
    double pos_arr[2][N];

    for(int i=0; i<N; i++){
        pos_arr[0][i] = 0;
        pos_arr[1][i] = 0;
    }

    double D=1, dt = 0.1, tmax = 100;
    double Exy = 0;
    double sigma = sqrt(2*D*dt);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(Exy, sigma);

    std::fstream dupa;
    dupa.open("dupa.dat", std::fstream::out);

    double tcondition = 0.1;
    for(double t=0; t<tmax; t+=dt){
        for(int i=0; i<N; i++){
            double dx = distribution(generator);
            double dy = distribution(generator);
            pos_arr[0][i] += dx;
            pos_arr[1][i] += dy;
        }
        if(t==tcondition){
            for(int i=0; i<N; i++){
                dupa<<pos_arr[0][i]<<" "<<pos_arr[1][i]<<"\n";
            }
            tcondition *= 10;
        }
    }
    dupa.close();

    return 0;

}