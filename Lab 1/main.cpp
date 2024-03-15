#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <time.h>

double rand_num_gen(){
    return (rand() / (double)RAND_MAX ); // pobiera ziarno w odwołaniu do funkcji void srand()
}

int bernoulli(double p){
    int X;
    if (rand_num_gen() <= p){
        X = 1;
    } else{
        X = 0;
    }
    
    return X;
}

int main(){
    srand((unsigned) time(NULL));

    double iterations = pow(10, 7);
    float propabilities[3] = {0.1, 0.5, 0.7};
    // std::cout<<propabilities[0]<<"  "<<propabilities[1]<<"  "<<propabilities[2]<<"  ";
    // for(int i=0; i<50; i++){std::cout<<bernoulli(propabilities[1])<<"\n";}
    std::fstream output;
    output.open("output.dat", std::fstream::out);

    for(int prop_i = 0; prop_i < 3; prop_i++){
        output<<propabilities[prop_i]<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<" "<<0<<"\n";

        double E_X_theo = propabilities[prop_i], var_theo = propabilities[prop_i]*(1-propabilities[prop_i]);

        int condition = 100;
        double X_sum = 0, X_sum_pow = 0;
        
        for(int i=0; i<iterations; i++){
            X_sum += bernoulli(propabilities[prop_i]);
            X_sum_pow += pow(bernoulli(propabilities[prop_i]),2);
            
            if(i+1 == condition){
                double N, E_X, E_X2, var;
                N = i+1;
                E_X = X_sum/(i+1);
                E_X2 = X_sum_pow/(i+1);
                var = E_X2 - pow(E_X,2);

                output<<N<<" "<<E_X<<" "<<E_X2<<" "<<var<<" "<<abs(E_X_theo-E_X)/E_X_theo<<" "<<abs(var_theo-var)/var_theo<<"\n";

                condition = condition * 10;
            }

        }
    }

    output.close();
}