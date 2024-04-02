#include "functions.h"

int main(){

    std::fstream log;
    log.open("log.dat", std::fstream::out);
    
    double ranges[3][3] = {
        {-3, 3}, // całka C1
        {0, 10}, // całka C2
        {0, 1} // całka C3
    };    
    double a, b;

    double (*func[3])(double) = {&fun_c_1, &fun_c_2, &fun_c_3};

    std::vector<double> con;

    int M = 10;

    double Nm_sum[6];

    double **Nm_arr = new double*[6];
    for(int i=0; i<6; i++){
        Nm_arr[i] = new double[M];
    } 
    
    // Metoda podstawowa

    std::fstream podst;
    podst.open("podst.dat", std::fstream::out);

    for(int c = 0; c <3; c++){
        a = ranges[c][0];
        b = ranges[c][1];

        basic_func(a, b, 2, 5, 1, func[c], &podst, &log, con);
    }
    podst.close();

    // Metoda losowania systematycznego (warstwowe nieoptymalne) 

    std::fstream sys;
    sys.open("sys.dat", std::fstream::out);

    for(int c = 0; c <3; c++){
        a = ranges[c][0];
        b = ranges[c][1];

        basic_func(a, b, 2, 5, M, func[c], &sys, &log, con);

        for(int j=0; j<M; j++){
            Nm_arr[2*c][j] = con[j];
                Nm_sum[2*c] += con[j];
            Nm_arr[2*c+1][j] = con[M+j];
                Nm_sum[2*c+1] += con[M+j];
        }
    }
    sys.close();

    // Metoda losowania warstwowego (optymalnego)

    for (int i=0; i<6; i++){
        for(int j=0; j<M; j++){
            Nm_arr[i][j] = Nm_arr[i][j]/Nm_sum[i];
            // std::cout<<Nm_arr[i][j]<<" ";
        }
        // std::cout<<"\n";
    }
    
    std::fstream opt;
    opt.open("opt.dat", std::fstream::out);

    for(int c = 0; c <3; c++){
        a = ranges[c][0];
        b = ranges[c][1];

        opt_func(a, b, 2, 5, M, Nm_arr, c, func[c], &opt, &log);
    }
    opt.close();

    log.close();

    return 0;
}