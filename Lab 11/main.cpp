#include "functions.h"
#include "photon_diffusion.h"


int main(){

    srand( time( NULL ) );

    time_t start, end;
    time(&start);

    for(int i=1; i<=8; i++){
    
        PHOTON_DIFFUSION_2D ob;
            ob.xmax = 0.2;
            ob.x_source = 0.1;
            ob.dx_source = 0.0;
            ob.x_detect = 0.15;
            ob.dx_detect = 0.01;
            ob.nx = 100;
            ob.ny = 100;
            ob.rx0 = 0.0;
            ob.ry0 = 1.0;
            ob.nlayers = 3;
            ob.layers_data[1][0] = 1.; // absorption
            ob.layers_data[1][1] = 10.; // scattering
            ob.layers_data[1][2] = 0.02; // width
            ob.layers_data[1][3] = 0.75; // g_anizo
            ob.layers_data[1][4] = 1.3; // n_refraction
            ob.layers_data[2][0] = 1.; // absorption
            ob.layers_data[2][1] = 190.; // scattering
            ob.layers_data[2][2] = 0.02; // width
            ob.layers_data[2][3] = 0.075; // g_anizo
            ob.layers_data[2][4] = 1.0; // n_refraction
            ob.layers_data[3][0] = 10.; // absorption
            ob.layers_data[3][1] = 90.; // scattering
            ob.layers_data[3][2] = 0.02; // width
            ob.layers_data[3][3] = 0.95; // g_anizo
            ob.layers_data[3][4] = 1.0; // n_refraction
        
        switch (i)
        {
        case 1:
            ob.rx0 = 0.8;
            ob.ry0 = 0.6;
            ob.layers_data[2][4] = 1.5;
            break;
        case 2:
            ob.rx0 = 0.8;
            ob.ry0 = 0.6;
            ob.layers_data[2][4] = 2.5;
            break;
        case 3:
            ob.rx0 = 0.8;
            ob.ry0 = 0.6;
            ob.layers_data[1][4] = 1.0;
            ob.layers_data[2][4] = 1.5;
            break;
        case 4:
            ob.rx0 = 0.8;
            ob.ry0 = 0.6;
            ob.layers_data[1][4] = 1.0;
            ob.layers_data[2][4] = 1.5;
            ob.layers_data[2][1] = 10;
            break;
        case 5:
            ob.rx0 = 0.0;
            ob.ry0 = 1.0;
            ob.layers_data[1][0] = 1.; // absorption
            ob.layers_data[1][1] = 10.; // scattering
            ob.layers_data[1][2] = 0.02; // width
            ob.layers_data[1][3] = 0.75; // g_anizo
            ob.layers_data[1][4] = 1.3; // n_refraction
            ob.layers_data[2][0] = 1.; // absorption
            ob.layers_data[2][1] = 190.; // scattering
            ob.layers_data[2][2] = 0.02; // width
            ob.layers_data[2][3] = 0.075; // g_anizo
            ob.layers_data[2][4] = 1.0; // n_refraction
            ob.layers_data[3][0] = 10.; // absorption
            ob.layers_data[3][1] = 90.; // scattering
            ob.layers_data[3][2] = 0.02; // width
            ob.layers_data[3][3] = 0.95; // g_anizo
            ob.layers_data[3][4] = 1.0; // n_refraction
            break;
        case 6:
            ob.layers_data[1][4] = 1.0; // n_refraction
            ob.layers_data[2][0] = 10.; // absorption
            ob.layers_data[2][1] = 210.; // scattering
            ob.layers_data[2][4] = 1.5; // n_refraction
            break;
        case 7:
            ob.layers_data[1][4] = 1.0; // n_refraction
            ob.layers_data[2][0] = 1.; // absorption
            ob.layers_data[2][1] = 210.; // scattering
            ob.layers_data[2][4] = 1.5; // n_refraction
            break;
        case 8:
            ob.layers_data[1][4] = 1.0; // n_refraction
            ob.layers_data[2][0] = 10.; // absorption
            ob.layers_data[2][1] = 210.; // scattering
            ob.layers_data[2][4] = 1.5; // n_refraction
            ob.layers_data[2][3] = 0.75; // g_anizo
            break;

        default:
            break;
        }
    
        ob.init();
        
        int N = 200000; // liczba wiązek fotonowych
        ob.write_all_paths = 0;
        ob.write_source_detection_paths = 0;

        for (int k=0; k<N; k++){
            ob.single_path();
        }

        std::string file_name = "absorption_sym_" + to_string(i) + ".dat";
        save_data(file_name, ob.absorption);

        std::vector<std::vector<double>> trans;
        trans.push_back(ob.transmittance);
        file_name = "transmittance_sym_" + to_string(i) + ".dat";
        save_data(file_name, trans);
        trans.clear();

        std::vector<std::vector<double>> refl;
        refl.push_back(ob.reflectance);
        file_name = "reflectance_sym_" + to_string(i) + ".dat";
        save_data(file_name, refl);
        refl.clear();
    }

    time(&end);
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(5); 
    std::cout << " sec \n"; 

    return 0;
}