#include "functions.h"


double uni_rand(){
    return (rand()/(double) RAND_MAX);
}

void save_data(std::string outFileName, std::vector<std::vector<double>> data){
    std::fstream outFile;
    outFile.open(outFileName, std::fstream::out);
    if(outFile.is_open()){
        for(int i=0; i < data.size(); i++){
            for(int j=0; j < data[i].size(); j++){
                outFile << data[i][j] << " ";
            }
            outFile << "\n";
        }
    }
    outFile.close();
}

std::vector<std::vector<double>> load_atoms(std::string inFileName){
    
    std::vector<std::vector<double>> xyz;

    std::ifstream inFile;
    inFile.open(inFileName.c_str());

    if(inFile.is_open()){
        while(!inFile.eof()){
            std::vector<double> row;
            row.resize(3,0);
            inFile >> row[0] >> row[1] >> row[2];
            xyz.push_back(row);
            // std::cout<<xyz[i][0]<<" "<<xyz[i][1]<<" "<<xyz[i][2]<<"\n";
        }
        inFile.close();
    }

    return xyz;
}

double rho_q(double rho_max, double sigma_p, int i, int j, int n_x, int n_y, double Delta){
    
    double r_x = i - Delta*n_x/2.0;
    double r_y = j - Delta*n_y/2.0;
    double r = pow(r_x, 2) + pow(r_y, 2);
    double rho = rho_max * exp(-1.0*r/(2*pow(sigma_p,2)));

    return rho;   
}

void relaksacja(const int n_x, const int n_y, double Delta, double epsilon, double omega, 
                double tol, double it_max, double rho_max, double sigma_p, double V_L,
                double V_T, double V_B){
    
    double F_old = 0;
    double F_new = 0;

    std::vector<std::vector<double>> V(n_x + 1, std::vector<double>(n_y + 1));
    std::vector<std::vector<double>> rho(n_x + 1, std::vector<double>(n_y + 1));
    
    for(int i=0; i <= n_y; i++){
        for(int j = 0; j <= n_x; j++){
            V[i][j] = 0;
            rho[i][j] = rho_q(rho_max, sigma_p, i, j, n_x, n_y, Delta);
            // std::cout<<"V = "<<V[i][j]<<"   rho = "<<rho[i][j]<<"\n";
        }
    }
    
    double x_max = n_x*Delta;
    double y_max = n_y*Delta;

    for(int j=0; j <= n_x; j++){
        V[j][0] = V_L * sin((M_PI*j)/(y_max));
    }
    for(int i=0; i <= n_y; i++){
        V[0][i] = V_B * sin((M_PI*i)/(x_max));
        V[n_y][i] = V_T * sin((M_PI*i)/(x_max));        
    }

    for(int it = 1; it < it_max; it++){
        for(int i=1; i<n_x; i++){
            for(int j=1; j<n_y; j++){
                double V_around = V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1];
                V[i][j] = (1-omega)*V[i][j] + omega/4.0 * (V_around + pow(Delta,2)/epsilon*rho[i][j]);
            }
        }

        for(int j=1; j<n_y; j++){
            V[j][n_x] = V[j][n_x-1];
        }

        F_old = F_new;
        F_new = 0;
        for(int i=1; i<n_y; i++){
            for(int j=1; j<n_x; j++){
                double E_y = (V[i+1][j]-V[i-1][j])/(2*Delta);
                double E_x = (V[i][j+1]-V[i][j-1])/(2*Delta);
                F_new += (pow(E_x,2) + pow(E_y,2))/2.0 - rho[i][j]*V[i][j];
            }
        }

        if(abs(1 - F_old/F_new) < tol){
            break;
        }
    }
    save_data("V_relaksacja.dat", V);
}

void MC_V(){
    
}