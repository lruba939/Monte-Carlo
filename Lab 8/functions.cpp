#include "functions.h"


double uni_rand(){
    return (rand()*0.99999/(double) RAND_MAX);
}

double cut_function(double r, double R1, double R2){
    double f_cut;
    if(r <= R1){
        f_cut = 1;
    } else if(r > R1 && r <= R2){
        double frac = (double) (r - R1)/(R2 - R1);
        f_cut = 1.0/2*(1 + cos( M_PI * frac ));
    } else{
        f_cut = 0;
    }
    return f_cut;
}

double Brenner(){

double R0, R1, R2, De, S, lambda, delta, a0, c0, d0;

R0 = 1.315;
R1 = 1.70;
R2 = 2.00;
De = 6.325;
S = 1.29;
lambda = 1.5;
delta = 0.80469;
a0 = 0.011304;
c0 = 19;
d0 = 2.5;


std::string inFileName = "data.dat";
int arrSize = 60;
double xyz[arrSize][3];

std::ifstream inFile;
inFile.open(inFileName.c_str());

    if(inFile.is_open()){
        for(int i=0; i<arrSize; i++){
            inFile >> xyz[i][0] >> xyz[i][1] >> xyz[i][2];
            // std::cout<<xyz[i][0]<<" "<<xyz[i][1]<<" "<<xyz[i][2]<<"\n";
        }
        inFile.close();
    }

    double Vi = 0;

    for(int i=0; i<arrSize; i++){
        for (int j = i+1; j < arrSize; j++)
        {
            // odleglosc miedzy atomami ij
            double xij, yij, zij;
            xij = xyz[i][0] - xyz[j][0];
            yij = xyz[i][1] - xyz[j][1];
            zij = xyz[i][2] - xyz[j][2];
            double rij[3] = {xij, yij, zij};
            double rij_len = sqrt( pow(xij,2) + pow(yij,2) + pow(zij,2) );

            // funkcja odcięcia ij
            double f_cut_ij = cut_function(rij_len, R1, R2);

            // potencjal odpychania ij
            double VR;
                VR = De/(S-1) * exp(-sqrt(2.*S)*lambda*(rij_len - R0));

            // potencjal przyciagania ij
            double VA;
                VA = De*S/(S-1) * exp(-sqrt(2./S)*lambda*(rij_len - R0));

            // czynnik skalujacy potencjal przyciagania

            double Bij = 0;
            double ksi = 0;
            double xik, yik, zik;
            for(int k=0; k<arrSize; k++){
                if(k != i){
                    xik = xyz[i][0] - xyz[k][0];
                    yik = xyz[i][1] - xyz[k][1];
                    zik = xyz[i][2] - xyz[k][2];
                    double rik[3] = {xik, yik, zik};
                    double rik_len = sqrt( pow(rik[0],2) + pow(rik[1],2) + pow(rik[2],2) );

                    double s_prod_ik = rij[0]*rik[0] + rij[1]*rik[1] + rij[2]*rik[2];
                    double len_ijk = rij_len * rik_len;
                    double kosinus = (s_prod_ik/len_ijk);
                    double g_function = a0*( 1 + pow(c0,2)/pow(d0,2) - pow(c0,2)/( pow(d0,2) + pow(1 + kosinus, 2) ) );
                    double f_cut_ik = cut_function(rik_len, R1, R2);
                    ksi += f_cut_ik*g_function;
                }
            }
            Bij = pow(1 + ksi, -delta);

            double Bji = 0;
            ksi = 0;
            double xjk, yjk, zjk;
            for(int k=0; k<arrSize; k++){
                if(k != j){
                    xjk = xyz[j][0] - xyz[k][0];
                    yjk = xyz[j][1] - xyz[k][1];
                    zjk = xyz[j][2] - xyz[k][2];
                    double rjk[3] = {xjk, yjk, zjk};
                    double rjk_len = sqrt( pow(rjk[0],2) + pow(rjk[1],2) + pow(rjk[2],2) );

                    double s_prod_jk = (double) rij[0]*rjk[0] + rij[1]*rjk[1] + rij[2]*rjk[2];
                    double len_ijk = rij_len * rjk_len;
                    double kosinus = (s_prod_jk/len_ijk);
                    double g_function = a0*( 1 + pow(c0,2)/pow(d0,2) - pow(c0,2)/( pow(d0,2) + pow(1 + kosinus, 2) ) );
                    double f_cut_jk = cut_function(rjk_len, R1, R2);
                    ksi += f_cut_jk*g_function;
                }
            }
            Bji = pow(1 + ksi, -delta);

            double B = (double) (Bij + Bji)/2;

            Vi += (double) f_cut_ij*( VR - B*VA );

        }
        
    }
    double V_tot = Vi/2.0;

    return V_tot;
}