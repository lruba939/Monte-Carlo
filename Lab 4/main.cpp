#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <math.h>


double normal_rand(){
    return (rand()/(double) RAND_MAX);
}

void kolo(double& X, double& Y, double R_mat, double x_r, double y_r){
    double U1, U2, U3;
    
    U1 = normal_rand();
    U2 = normal_rand();
    U3 = normal_rand();

    X = sqrt(-2.*log(1-U1))*cos(2*M_PI*U2);
    Y = sqrt(-2.*log(1-U1))*sin(2*M_PI*U2);

    double r = sqrt(X*X + Y*Y);

    X = sqrt(U3) * X/r * R_mat + x_r;
    Y = sqrt(U3) * Y/r * R_mat + y_r;
}

int czesc_wspolna(double X, double Y, double x_a, double y_a, double x_b, double y_b, double R_A, double R_B){
    if( ( sqrt(pow((X-x_a),2) + pow((Y-y_a),2)) <= R_A) && ( sqrt(pow((X-x_b),2) + pow((Y-y_b),2)) <= R_B) ){
        return 1;
    }else{
        return 0;
    }
}


int main(){

    int N = pow(10,4);

    std::fstream KA; // dane dla koła o parametrach A
    KA.open("KA.dat", std::fstream::out);
    std::fstream KB; // dane dla koła o parametrach B
    KB.open("KB.dat", std::fstream::out);

    double R_A = 1, x_a = 0, y_a = 0;
    double R_B = 1, x_b = 0, y_b = 0;

    R_A = 1;
    R_B = sqrt(2)*R_A;
    x_a = R_A + R_B;

    for(int i=0; i<N; i++){
        double X = 0, Y = 0;
        kolo(X, Y, R_A, x_a, y_a); // punkty dla kola A
        KA<<X<<" "<<Y<<"\n";

        X = 0, Y = 0;
        kolo(X, Y, R_B, x_b, y_b); // punkty dla kola B
        KB<<X<<" "<<Y<<"\n";
    }

    KA.close();
    KB.close();

    // część odpowiedzialna za całkowanie

    R_A = 1; x_a = 0; y_a = 0;
    R_B = sqrt(2)*R_A; x_b = 0; y_b = 0;


    std::fstream data;
    // podpunkt a)
    // double R = R_A; x_a = R_B + 0.5*R_A;
    //     data.open("3a.dat", std::fstream::out);

    // podupunkt b)
    // double R = R_A; x_a = 0;
    //     data.open("3b.dat", std::fstream::out);

    // podupnkt c)
    // double R = R_B; x_a = R_B + 0.5*R_A;
    //     data.open("3c.dat", std::fstream::out);

    // podupunkt d)
    double R = R_B; x_a = 0;
        data.open("3d.dat", std::fstream::out);

    for(int p=2; p<=6; p++){
        double suma = 0;
        N = pow(10,p);

        for (int i = 0; i < N; i++)
        {
            double X = 0, Y = 0;
            kolo(X, Y, R, x_a, y_a);
            suma += czesc_wspolna(X, Y, x_a, y_a, x_b, y_b, R_A, R_B);
        }

        std::cout<<suma<<"\n";
        
        double mu1 = M_PI*R*R*suma/N;
        double mu2 = M_PI*R*R*mu1;
        double var = (mu2 - mu1*mu1)/N;
        double std = sqrt(var);

        data<<p<<" "<<mu1<<" "<<mu2<<" "<<var<<" "<<std<<"\n";
    }

    data.close();

    return 0;
}