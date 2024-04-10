#include "functions.h"

double uni_rand(){
    return (rand()*0.99999/(double) RAND_MAX);
}

void normal_rand(double& X, double& Y, double std){
    double U1, U2;
    
    U1 = uni_rand();
    U2 = uni_rand();

    X = sqrt(-2.*log(1-U1))*cos(2*M_PI*U2)*std;
    Y = sqrt(-2.*log(1-U1))*sin(2*M_PI*U2)*std;
}

void wiener(int N, double tmax, double dt, double D, std::fstream *pos, std::fstream *stat){
    double pos_arr[2][N];
  
    double sigma = sqrt(2*D*dt);

    for(int i=0; i<N; i++){
        pos_arr[0][i] = 0;
        pos_arr[1][i] = 0;
    }

    double dx, dy;
    double tcondition = 0.1;
    for(double t=0; t<=tmax; t+=dt){
        std::cout<<"Particles = "<<N<<" time = "<<t<<"\n";

        std::vector<double> stat_values_x = {0,0,0}; // {mean, second moment, D_XX}
        std::vector<double> stat_values_y = {0,0,0}; // {mean, second moment, D_YY}
        std::vector<double> stat_values_xy = {0,0,0}; // {mean, second moment, D_XY}

        for(int i=0; i<N; i++){
            normal_rand(dx, dy, sigma);
            pos_arr[0][i] += dx;
            pos_arr[1][i] += dy;

            stat_values_x[0] = stat_values_x[0] + pos_arr[0][i]/N; // X mean
            stat_values_x[1] += pow(pos_arr[0][i],2)/N; // X second moment
            
            stat_values_y[0] += pos_arr[1][i]/N; // Y mean
            stat_values_y[1] += pow(pos_arr[1][i],2)/N; // Y second moment

            stat_values_xy[1] += pos_arr[0][i]*pos_arr[1][i]/N; // XY second moment
        }
        stat_values_xy[0] = stat_values_x[0]*stat_values_y[0]; // XY mean

        if(t!=0){
            stat_values_x[2] = (stat_values_x[1] - pow(stat_values_x[0],2))/(2.*t); // D_XX
            stat_values_y[2] = (stat_values_y[1] - pow(stat_values_y[0],2))/(2.*t); // D_YY
            stat_values_xy[2] = (stat_values_xy[1] - stat_values_xy[0])/(2.*t); // D_XY
        }else{
            stat_values_x[2] = NAN;
            stat_values_y[2] = NAN;
            stat_values_xy[2] = NAN;
        }

        *stat<<N<<" "<<t<<" "<<stat_values_x[0]<<" "<<stat_values_y[0]<<" "<<stat_values_xy[0]<<" "<<stat_values_x[2]<<" "<<stat_values_y[2]<<" "<<stat_values_xy[2]<<"\n";

        if(int(std::ceil(t*10))==int(std::ceil(tcondition*10))){
            std::cout<<"Done!\n";
            for(int i=0; i<N; i++){
                    *pos<<pos_arr[0][i]<<" "<<pos_arr[1][i]<<" "<<t<<"\n";
                }
                tcondition *= 10;
        }
    }
}


void dyf_and_abs(int N, double tmax, double dt, double xs, double ys, double D, double omega, double xr, double yr, double Rr, double xa, double ya, double Ra, std::fstream *pos, std::fstream *stat){
    double pos_arr[3][N];
  
    double sigma = sqrt(2*D*dt);

    for(int i=0; i<N; i++){
        pos_arr[0][i] = 0; // x
        pos_arr[1][i] = 0; // y
        pos_arr[2][i] = 0; // aktywnosc czastek 
    }

    double tcondition = 0.1;
    for(double t=0; t<=tmax; t+=dt){
        int n = 0, nnew = 0;

        double dn = omega * t;

        std::vector<double> stat_values_x = {0,0,0}; // {mean, second moment, D_XX}
        std::vector<double> stat_values_y = {0,0,0}; // {mean, second moment, D_YY}
        std::vector<double> stat_values_xy = {0,0,0}; // {mean, second moment, D_XY}

        for(int i=1; i<=N; i++){
            
            if(pos_arr[2][i]==0 && nnew<dn){
                pos_arr[2][i] = 1;
                pos_arr[0][i] = xs;
                pos_arr[1][i] = ys;
                nnew++;
            }

            double x1, y1, x2, y2, dx, dy;
            if(int(pos_arr[2][i])==1){
                x1 = pos_arr[0][i];
                y1 = pos_arr[1][i];

                normal_rand(dx, dy, sigma);
                x2 = x1 + dx;
                y2 = y1 + dy;

                double length = 0;
                int theta = pos_arr[2][i];
                do{
                    particle_translation(x1, y1, x2, y2, xr, yr, Rr, xa, ya, Ra, theta, length);
                }while(length > pow(10, -6));
                pos_arr[2][i] = theta;

                pos_arr[0][i] = x2;
                pos_arr[1][i] = y2;

                if(pos_arr[2][i]){
                    n++;
                }
            }
            stat_values_x[0] = stat_values_x[0] + pos_arr[0][i]/N; // X mean
            stat_values_x[1] += pow(pos_arr[0][i],2)/N; // X second moment
            
            stat_values_y[0] += pos_arr[1][i]/N; // Y mean
            stat_values_y[1] += pow(pos_arr[1][i],2)/N; // Y second moment

            stat_values_xy[1] += pos_arr[0][i]*pos_arr[1][i]/N; // XY second moment
        }
        stat_values_xy[0] = stat_values_x[0]*stat_values_y[0]; // XY mean

        if(t!=0){
            stat_values_x[2] = (stat_values_x[1] - pow(stat_values_x[0],2))/(2.*t); // D_XX
            stat_values_y[2] = (stat_values_y[1] - pow(stat_values_y[0],2))/(2.*t); // D_YY
            stat_values_xy[2] = (stat_values_xy[1] - stat_values_xy[0])/(2.*t); // D_XY
        }else{
            stat_values_x[2] = NAN;
            stat_values_y[2] = NAN;
            stat_values_xy[2] = NAN;
        }
        *stat<<N<<" "<<t<<" "<<stat_values_x[0]<<" "<<stat_values_y[0]<<" "<<stat_values_xy[0]<<" "<<stat_values_x[2]<<" "<<stat_values_y[2]<<" "<<stat_values_xy[2]<<" "<<n<<"\n";

        std::cout<<t<<" "<<t/tmax*100<<"%\n";
        if(int(std::ceil(t*10))==int(std::ceil(tcondition*10))){
            std::cout<<"Done!\n";
            for(int i=0; i<N; i++){
                    *pos<<pos_arr[0][i]<<" "<<pos_arr[1][i]<<" "<<pos_arr[2][i]<<" "<<t<<"\n";
                }
            tcondition *= 10;
        }
    }
}