#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <algorithm>


double zlozlony(double g1){
    double U1, U2, X;

    U1 = rand()/(double)RAND_MAX;
    U2 = rand()/(double)RAND_MAX;

    if(U1 <= g1){
        X = U2;
    } else{
        X = sqrt( 1 - sqrt(1-U2) );
    }

    return X;
}

double markowa(double X, double Delta){
    double U1, U2;

    U1 = rand()/(double)RAND_MAX;
    U2 = rand()/(double)RAND_MAX;

    double X_i = X, X_ii;
    double x_new = X_i + (2*U1-1)*Delta;

    double p_acc = std::min((double) 1, (double) (1 + x_new - pow(x_new,3))/(1 + X_i - pow(X_i,3)));

    if(x_new <= 1 && x_new >= 0 && U2 <= p_acc){
        X_ii = x_new;
    }else{
        X_ii = X_i;
    }
    return X_ii;
}

double eliminacji(double G){
    double U1, U2, G2;
    do{
        U1 = rand()/(double)RAND_MAX;
        U2 = rand()/(double)RAND_MAX;
        G2 = G*U2;
    }while ((G2 > ( (double) 4/5*(1 + U1 - pow(U1,3)) ) ));
    
    return U1;
}

int main(){
    // srand(time ( NULL ));

    int N = pow(10, 6);

    double g1 = (double) 4/5;

    // Złożony
    std::fstream outputZlozony;
        outputZlozony.open("output_zlozony.dat", std::fstream::out);

    for(int i = 0; i < N; i++){
        outputZlozony<<i+1<<" "<<zlozlony(g1)<<"\n";

    }   outputZlozony.close();

    // Markowa
    std::fstream outputMarkowa1;
        outputMarkowa1.open("output_Markowa1.dat", std::fstream::out);
    std::fstream outputMarkowa2;
        outputMarkowa2.open("output_Markowa2.dat", std::fstream::out);

    double X_i1 = rand()/(double)RAND_MAX;
    double X_i2 = rand()/(double)RAND_MAX;
    for(int i = 0; i < N; i++){
        double X1 = markowa(X_i1, 0.5);
        outputMarkowa1<<i+1<<" "<<X1<<"\n";
        double X2 = markowa(X_i2, 0.05);
        outputMarkowa2<<i+1<<" "<<X2<<"\n";
        X_i1 = X1;
        X_i2 = X2;

    }   outputMarkowa1.close(); outputMarkowa2.close();

        double G = 1.15;

    // Eliminacji
    std::fstream outputEliminacji;
        outputEliminacji.open("output_eliminacji.dat", std::fstream::out);

    for(int i = 0; i < N; i++){
        outputEliminacji<<i+1<<" "<<eliminacji(G)<<"\n";

    }   outputEliminacji.close();

    return 0;
}