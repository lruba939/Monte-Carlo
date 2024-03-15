#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <math.h>


double normal_rand(){
    return (rand()/(double) RAND_MAX);
}

void boxa_mullera(double& X, double& Y){
    double U1, U2;
    
    U1 = normal_rand();
    U2 = normal_rand();

    X = sqrt(-2.*log(1-U1))*cos(2*M_PI*U2);
    Y = sqrt(-2.*log(1-U1))*sin(2*M_PI*U2);
}


int main(){

    std::fstream BM;
    BM.open("boxa-mullera.dat", std::fstream::out);

    std::fstream Kolo;
    Kolo.open("kolo.dat", std::fstream::out);

    // definicja wartości dla transformacji afinicznej
    double b1 = 1, b2 = 0.2;
    double alpha = M_PI/4;
    double kosinus = cos(alpha), sinus = sin(alpha);
    double sum_elipX = 0, sum_elipY = 0, sum2_elipX = 0, sum2_elipY = 0, sum_elipXY = 0;
    double BMsum_elipX = 0, BMsum_elipY = 0, BMsum2_elipX = 0, BMsum2_elipY = 0, BMsum_elipXY = 0;

    std::fstream elipsa;
    elipsa.open("elipsa.dat", std::fstream::out);
    std::fstream mat_cov;
    mat_cov.open("mat_cov.dat", std::fstream::out);
    std::fstream elipsa_norm;
    elipsa_norm.open("elipsa_norm.dat", std::fstream::out);
    std::fstream mat_cov_norm;
    mat_cov_norm.open("mat_cov_norm.dat", std::fstream::out);


    int N = pow(10,4);
    double X = 0, Y = 0;

    for(int i = 0; i <= N; i++){

        // fragment generujący wartości dla modelu Boxa-Mullera
        boxa_mullera(X, Y);
            BM<<X<<" "<<Y<<"\n";

        // fragment używający poprzednich wartości X i Y do stworzenia rozkładu koła jednostkowego o rozkładzie normalnym
        double dX = X / ( sqrt( X*X + Y*Y ) );
        double dY = Y / ( sqrt( X*X + Y*Y ) );
        double R = sqrt(normal_rand());
        double ddX = dX*R, ddY = dY*R;
            Kolo<<ddX<<" "<<ddY<<"\n";

        // fragment dla utworzenia elipsy na podstawie wyników koła

        double elipX, elipY;
        elipX = ddX*b1*kosinus - ddY*b2*sinus;
        elipY = ddX*b1*sinus + ddY*b2*kosinus;
            elipsa<<elipX<<" "<<elipY<<"\n";
        
        sum_elipX += elipX;
            sum2_elipX += elipX*elipX;
        sum_elipY += elipY;
            sum2_elipY += elipY*elipY;
        sum_elipXY += elipX*elipY;

        // fragment dla utworzenia elipsy na podstawie wyników Boxa-Mullera

        double BMelipX, BMelipY;
        BMelipX = X*b1*kosinus - Y*b2*sinus;
        BMelipY = X*b1*sinus + Y*b2*kosinus;
            elipsa_norm<<BMelipX<<" "<<BMelipY<<"\n";
        
        BMsum_elipX += BMelipX;
            BMsum2_elipX += BMelipX*BMelipX;
        BMsum_elipY += BMelipY;
            BMsum2_elipY += BMelipY*BMelipY;
        BMsum_elipXY += BMelipX*BMelipY;
    }

    double sigma2_elipX = sum2_elipX/N - pow(sum_elipX/N,2);
    double sigma2_elipY = sum2_elipY/N - pow(sum_elipY/N,2);
    double sigma_elipXY = sum_elipXY/N - sum_elipX/N*sum_elipY/N;
    double r_xy = sigma_elipXY / ( sqrt( sigma2_elipX*sigma2_elipY ) );
        mat_cov<<sigma2_elipX<<" "<<sigma_elipXY<<"\n";
        mat_cov<<sigma_elipXY<<" "<<sigma2_elipY<<"\n";
        mat_cov<<r_xy;

    double BMsigma2_elipX = BMsum2_elipX/N - pow(BMsum_elipX/N,2);
    double BMsigma2_elipY = BMsum2_elipY/N - pow(BMsum_elipY/N,2);
    double BMsigma_elipXY = BMsum_elipXY/N - BMsum_elipX/N*BMsum_elipY/N;
    double BMr_xy = BMsigma_elipXY / ( sqrt( BMsigma2_elipX*BMsigma2_elipY ) );
        mat_cov_norm<<BMsigma2_elipX<<" "<<BMsigma_elipXY<<"\n";
        mat_cov_norm<<BMsigma_elipXY<<" "<<BMsigma2_elipY<<"\n";
        mat_cov_norm<<BMr_xy;


    BM.close();
    Kolo.close();
    elipsa.close();
    mat_cov.close();
    elipsa_norm.close();
    mat_cov_norm.close();

    return 0;
}