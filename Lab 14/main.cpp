#include "utils.h"
#include "vqmc.h"


int main(){

    // srand( time( NULL ) );

    time_t start, end;
    time(&start);


    double a = 1.0, c = 0.0, r = 0.1, dr = 0.1, N = pow(10,6);

    std::vector<std::vector<double>> histogram_data;
    histogram_data.push_back(histogram(N, dr, r, a, c));
    save_data("histogram.dat", histogram_data);

    double a_start = 0.3, a_stop = 1.2, da = 0.02,
            c_start = -0.7, c_stop = 0.3, dc = 0.02;

    std::vector<std::vector<double>> map_en, map_std;
    std::vector<double> row_en, row_std, con;

    r = 0.1;
    c = c_stop;
    while(c >= c_start){
        a = a_stop;
        while(a >= a_start){
            con = energy(N, dr, r, a, c);
            row_en.push_back(con[0]);
            row_std.push_back(sqrt(con[1]));
            a -= da;
        }
        map_en.push_back(row_en);
        map_std.push_back(row_std);
        row_en.clear();
        row_std.clear();
        c -= dc;
    }

    save_data("map_en.dat", map_en);
    save_data("map_std.dat", map_std);


    time(&end);
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(5); 
    std::cout << " sec \n"; 

    return 0;
}