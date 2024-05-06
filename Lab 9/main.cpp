#include "functions.h"


int main(){

    srand( time( NULL ) );

    time_t start, end;
    time(&start);

    int n_x = 30;
    int n_y = 30;
    double Delta = 0.1;
    double V_L = 1.0;
    double V_T = -1.0;
    double V_B = -1.0;
    double epsilon = 1.0;
    double x_max = Delta*n_x;
    double y_max = Delta*n_y;
    double rho_max = 1.0;
    double sigma_p = x_max/10.0;

    int it_max = pow(10,4);
    double tol = pow(10, -6);
    double omega = 1.8;

    relaksacja(n_x, n_y, Delta, epsilon, omega, tol, it_max, rho_max, sigma_p, V_L, V_T, V_B);

    //////////////////////////
    int N_chains_max = 300;
    double n_length_max = 300;
    int B_condition = 1;
    MC_V(n_x, n_y, Delta, N_chains_max, n_length_max, epsilon, rho_max, sigma_p, V_L, V_T, V_B, B_condition);

    time(&end);
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(5); 
    std::cout << " sec \n"; 

    return 0;
}