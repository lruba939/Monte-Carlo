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

std::vector<std::vector<double>> load_data(std::string inFileName){
    
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

double R_0_f(double L, double C){
    double R_0 = sqrt(L/C);

    return R_0;
}

double zeta_f(double R_0, double R_g){
    double zeta = R_0 / (R_0 + R_g);

    return zeta;
}

double Gamma_g_f(double R_0, double R_g){
    double Gamma_g = (R_g - R_0) / (R_g + R_0);

    return Gamma_g;
}

double Gamma_l_f(double R_0, double R_l){
    double Gamma_l = (R_l - R_0) / (R_l + R_0);

    return Gamma_l;
}

double c_f(double L, double C){
    double c = sqrt(1 / (L*C));

    return c;
}

double mu_f(double G, double C){
    double mu = G / C;

    return mu;
}

double lambda_f(double R, double L, double G, double C){
    double lambda = (R/L - G/C)/2.0;

    return lambda;
}

double V_g(double t, double s){
    double nu = 1 * pow(10,9);
    double t_0 = 7.5 * pow(10,-9);
    double sigma = 0.75 * pow(10,-9);

    double V_g = sin(2*M_PI*nu*t)*exp( (-1) * pow((t-t_0),2) / (2*pow(sigma,2)) );

    return V_g;
}

double MC_line(int npaths, double x_start, double t_start, double c, double lambda, double mu,
             double zeta, double Gamma_g, double Gamma_l, double l){

    double fxt = 0.;
    double bxt = 0.;

    for(int i=1; i<=2; i++){
        double suma = 0.;

        for(int n=0; n<npaths; n++){
            double x = x_start;
            double t = t_start;
            double eta = 1.0;
            int sign = pow(-1, i);

            while(t > 0){
                double U1 = uni_rand();
                double s = (-1)*(log(U1)) / (lambda + mu);
                
                if(sign == -1){
                    if((x - c*s) > 0){
                        eta = eta * lambda / (lambda + mu);
                    } else{
                        s = x / c;
                        suma += eta*zeta*V_g(t-s, s);
                        eta = eta*Gamma_g;
                    }
                    x = x-c*s;
                    t = t-s;
                }
                if(sign == 1){
                    if((x+c*s) < l){
                        eta = eta * lambda / (lambda + mu);
                    } else{
                        s = (l-x)/c;
                        eta = eta*Gamma_l;
                    }
                    x = x+c*s;
                    t = t-s;
                }
                sign = (-1)*sign;
            }
        }

        if(i==1){
            fxt = suma/npaths;
        }
        if(i==2){
            bxt = suma/npaths;
        }
    }

    double uxt = fxt + bxt;

    return uxt;
}