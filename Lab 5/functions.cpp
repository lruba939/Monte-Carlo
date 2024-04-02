#include "functions.h"

double normal_rand(){
    return (rand()/(double) RAND_MAX);
}

double rand_in_range(double a, double b){
    return a + (b - a)*normal_rand();
}

double fun_c_1(double x){
    return 1. + tanh(x);
}

double fun_c_2(double x){
    return 1./(1+pow(x,2));
}

double fun_c_3(double x){
    return pow(cos(M_PI*x),10);
}

void basic_func(double a, double b, int k_start, int k_stop, int M, double (*function_C)(double), std::fstream *infile, std::fstream *log, std::vector<double>& con){
    double xi = 0; double result = 0;
    double delta_x = (b - a)/M;
    double pm = 1./M;
    std::vector<double> std_m;
    
    for(int k = k_start; k <= k_stop; k++){
        int N = pow(10,k)*pm;
        std::vector<double> stat_values = {0,0,0,0,0}; // {mean, second moment, var, std, error}

        for(int m=0; m<M; m++){
            double x_start = a + delta_x*m;
            double x_stop = x_start + delta_x;

            for(int i=0; i<N; i++){
                xi = rand_in_range(x_start, x_stop);
                if(k==2 && b==3){
                    *log<<xi<<"\n";
                }
                result = (b - a)*function_C(xi);
                stat_values[0] += result/N;
                stat_values[1] += pow(result, 2)/N;
                stat_values[2] = (pow(result, 2)/N - pow(result/N, 2));
            }
            if(k==2 || k==3){
                std_m.push_back(sqrt(stat_values[2]/N));
            }
        }
        con = std_m;

        stat_values[0] = stat_values[0]*pm;
        stat_values[1] = stat_values[1]*pm;
        stat_values[2] = (stat_values[1] - pow(stat_values[0], 2))/N*pm*pm;
        stat_values[3] = sqrt(stat_values[2]);
        stat_values[4] = stat_values[3]/stat_values[0]*100;

        *infile<<k<<" "<<stat_values[0]<<" "<<stat_values[1]<<" "<<stat_values[2]<<" "<<stat_values[3]<<" "<<stat_values[4]<<"\n";
    }
    std::cout<<"Done. \n";
}


void opt_func(double a, double b, int k_start, int k_stop, int M, double **Nm_arr, int c, double (*function_C)(double), std::fstream *infile, std::fstream *log){
    double xi = 0; double result = 0;
    double delta_x = (b - a)/M;
    double pm = 1./M;
    int N = 1;
    
    for(int k = k_start; k <= k_stop; k++){
        std::vector<double> stat_values = {0,0,0,0,0}; // {mean, second moment, var, std, error}
        
        for(int m=0; m<M; m++){
            if(k==2){
                // std::cout<<Nm_arr[c*2][m]<<" \n";
                N = std::ceil( Nm_arr[c*2][m]*pow(10,k) );
            }else{
                // std::cout<<Nm_arr[c*2+1][m]<<" \n";
                N = std::ceil( Nm_arr[c*2+1][m]*pow(10,k) );
            }
            // std::cout<<N<<" ";
            double x_start = a + delta_x*m;
            double x_stop = x_start + delta_x;

            for(int i=0; i<N; i++){
                xi = rand_in_range(x_start, x_stop);
                if(k==2 && b==3){    
                    *log<<xi<<"\n";
                }
                result = (b - a)*function_C(xi);
                stat_values[0] += result/N;
                stat_values[1] += pow(result, 2)/N;
            }
        }
        stat_values[0] = stat_values[0]*pm;
        stat_values[1] = stat_values[1]*pm;
        stat_values[2] = (stat_values[1] - pow(stat_values[0], 2))/N*pm*pm;
        stat_values[3] = sqrt(stat_values[2]);
        stat_values[4] = stat_values[3]/stat_values[0]*100;

        *infile<<k<<" "<<stat_values[0]<<" "<<stat_values[1]<<" "<<stat_values[2]<<" "<<stat_values[3]<<" "<<stat_values[4]<<"\n";
    }
    std::cout<<"Done. \n";
}