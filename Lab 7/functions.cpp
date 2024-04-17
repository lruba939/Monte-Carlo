#include "functions.h"

double uni_rand(){
    return (rand()*0.99999/(double) RAND_MAX);
}


void Gillespie(double k1, double k2, double k3, double k4, double x1start, double x2start, double x3start, int tmax, int N, int Pmax, std::fstream *stat, std::fstream *hist){

    double h0[N], h1[N], h2[N];
    double dt = tmax/N;
    int ncount[N];

    for(int i=0; i<N; i++){
        h0[i] = 0;
        h1[i] = 0;
        h2[i] = 0;
        ncount[i] = 0;
    }

    for(int p=1; p<=Pmax; p++){
        double t = 0;
        double x1 = x1start;
        double x2 = x2start;
        double x3 = x3start;

        while(t<tmax){
            double G1 = k1;
            double G2 = k2;
            double G3 = k3*x1*x2;
            double G4 = k4*x3;

            double Gmax = G1+G2+G3+G4;

            double U1 = uni_rand();
            double U2 = uni_rand();

            double Dt = (-1)/Gmax*log(U1);

            if(U2*Gmax <= G1){
                x1 += 1;
            }else if(U2*Gmax <= G1+G2){
                x2 += 1;
            }else if(U2*Gmax <= G1+G2+G3){
                x1 -= 1;
                x2 -= 1;
                x3 += 1;
            }else{
                x3 -= 1;
            }

            *stat<<t+Dt<<" "<<x1<<" "<<x2<<" "<<x3<<"\n";

            int l = floor(t/dt);
            // std::cout<<h0[l]<<" ";
            h0[l] += x3;
            ncount[l] += 1;

            // std::cout<<"count: "<<ncount[l]<<"   val: "<<h0[l]<<"   x3: "<<x3<<"\n";

            t += Dt;

            // std::cout<<"P: "<<p<<"   Time: "<<t<<"   Perc: "<<t/tmax*100<<"\n";
        }

        for(int i=0; i<N; i++){
            double avg_x3 = h0[i]/ncount[i];
            h1[i] += avg_x3;
            h2[i] += pow(avg_x3,2);
        }

    }
    
    double hstd[N];
    // HIST
    for(int i=0; i<N; i++){
        double x3t_1 = h1[i]/Pmax;
        double x3t_2 = h2[i]/Pmax;
        hstd[i] = sqrt( (x3t_2 - pow(x3t_1, 2)) / Pmax );
        // std::cout<<i<<"\n";
        *hist<<i*dt+dt/2<<" "<<x3t_1<<" "<<hstd[i]<<"\n";
    }

}
