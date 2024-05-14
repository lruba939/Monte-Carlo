#include "functions.h"
#include "u_xt_exact.h"

int main(){

    srand( time( NULL ) );

    time_t start, end;
    time(&start);

    // Parametry symulacji
    double L = 0.25 * pow(10, -6);
    double C = 100. * pow(10, -12);
    double R = 12.5;
    double G = 0.5 * pow(10, -3);
    double l = 2.0;
    double R_l = 12.5;
    double R_g = 75.;

    double R_0 = R_0_f(L, C);
    double zeta = zeta_f(R_0, R_g);
    double Gamma_g = Gamma_g_f(R_0, R_g);
    double Gamma_l = Gamma_l_f(R_0, R_l);

    double c = c_f(L, C);
    double lambda = lambda_f(R, L, G, C);
    double mu = mu_f(G, C);

    std::vector<std::vector<double>> line;
    std::vector<double> row;
    row.resize(3,0);

    int npaths = pow(10,3);
    double t = 22. * pow(10,-9);
    double uxt = 0;
    double uxt_theo = 0;

    double delta = 0.01;
    double x = 0;
    for(int i=0; i<=l/delta; i++){
        uxt = MC_line(npaths, x, t, c, lambda, mu, zeta,
                      Gamma_g, Gamma_l, l);
        uxt_theo = u_xt_exact(x, t, 7.5*pow(10,-9), 1*pow(10,9), 0.75*pow(10,-9), R, G, L, C, R_g, R_l, l, 1000, 100);

        row[0] = x;
        row[1] = uxt;
        row[2] = uxt_theo;
        line.push_back(row);
        std::cout<<"x = "<<x<<"\n";
        x = x + delta;
    }

    save_data("test.dat", line);


    time(&end);
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(5); 
    std::cout << " sec \n"; 

    return 0;
}