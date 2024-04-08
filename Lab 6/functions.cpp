#include "functions.h"

double uni_rand(){
    return (rand()/(double) RAND_MAX);
}

void normal_rand(double& X, double& Y, double mean, double std){
    double U1, U2;
    
    U1 = uni_rand();
    U2 = uni_rand();

    X = mean + sqrt(-2.*log(1-U1))*cos(2*M_PI*U2)*std;
    Y = mean + sqrt(-2.*log(1-U1))*sin(2*M_PI*U2)*std;


}

void wiener(int N, double tmax, double dt, double Exy, double D, std::fstream *infile){
    double pos_arr[2][N];

    double sigma = sqrt(2*D*dt);

    for(int i=0; i<N; i++){
        pos_arr[0][i] = 0;
        pos_arr[1][i] = 0;
    }

    double dx, dy;
    double tcondition = 0.1;
    for(double t=0; t<tmax; t+=dt){
        for(int i=0; i<N; i++){
            normal_rand(dx, dy, Exy, sigma);
            pos_arr[0][i] += dx;
            pos_arr[1][i] += dy;
        }
        
    std::cout<<"Particles = "<<N<<" time = "<<t<<"\n";

    if(int(std::ceil(t*10))==int(std::ceil(tcondition*10))){
        std::cout<<"Done!\n";
        for(int i=0; i<N; i++){
                *infile<<pos_arr[0][i]<<" "<<pos_arr[1][i]<<" "<<t<<"\n";
            }
            tcondition *= 10;
        }
    }
}