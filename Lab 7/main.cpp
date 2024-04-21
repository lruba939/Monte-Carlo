#include "functions.h"



int main(){

    //std::srand(time(NULL));

    double k1, k2, k3, k4, x1start, x2start, x3start;
    int tmax, N, Pmax;

    k1 = 1.;
    k2 = 1.;
    k3 = 0.001;
    k4 = 0.01;

    x1start = 120; // for t=0
    x2start = 80; // for t=0
    x3start = 1; // for t=0

    tmax = 200;
    N = 50;
    Pmax = 5;

    std::fstream stat, hist;
    stat.open("stat.dat", std::fstream::out);
    hist.open("hist.dat", std::fstream::out);

        Gillespie(k1, k2, k3, k4, x1start, x2start, x3start, tmax, N, Pmax, &stat, &hist);

    stat.close();
    hist.close();


    Pmax = 100;

    stat.open("stat2.dat", std::fstream::out);
    hist.open("hist2.dat", std::fstream::out);

        Gillespie(k1, k2, k3, k4, x1start, x2start, x3start, tmax, N, Pmax, &stat, &hist);

    stat.close();
    hist.close();

    return 0;

}