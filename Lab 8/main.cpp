#include "functions_2.h"


int main(){

    srand( time( NULL ) );

    time_t start, end;
    time(&start);

    double V_tot = 0;

    // Liczenie wartości teoretycznych dla fullerenu
    std::vector<std::vector<double>> fullerene_file_data = load_atoms("data.dat");
    std::vector<std::vector<double>> fullerene;
    std::vector<double> atoms_row(6);
    for (int i = 0; i < fullerene_file_data.size(); i++){
        atoms_row[0] = 0;
        atoms_row[1] = 0;
        atoms_row[2] = 0;
        atoms_row[3] = fullerene_file_data[i][0];
        atoms_row[4] = fullerene_file_data[i][1];
        atoms_row[5] = fullerene_file_data[i][2];
        fullerene.push_back(atoms_row);
    }
    
    for(int i = 0; i < fullerene_file_data.size(); i++){
        V_tot += Brenner_i_energy(fullerene, i, false)/2.0;
    }

    std::cout<<"Theoretical fullerene total energy: "<<V_tot<<" eV \n";

    // Parametry podstawowe symulacji
    const int n = 60;
    const int it_max = pow(10,5);
    int M = 100;
    double w_r = pow(10.0,-4);
    double w_phi = pow(10.0,-2)*5;
    double w_theta = pow(10.0,-2)*5;
    double W_all = pow(10.0,-4);
    double r_start = 3.5;
    double beta_min = 1.0;
    double beta_max = 100.0;
    double p = 2.0;
    bool correction = false;

    std::ofstream energyInFile, posInFile, pcfInFile, avogadroInFile;

    // // Zadanie trzecie
    // energyInFile.open("atoms_energy_zad3.dat");
    // posInFile.open("atoms_pos_zad3.dat");
    // pcfInFile.open("pcf_zad3.dat");
    // avogadroInFile.open("atoms_avogadro_zad3.xyz");

    // SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
    //     beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

    // energyInFile.close();
    // posInFile.close();
    // pcfInFile.close();
    // avogadroInFile.close();

    // // Zadanie czwarte
    // correction = true;

    // energyInFile.open("atoms_energy_zad4.dat");
    // posInFile.open("atoms_pos_zad4.dat");
    // pcfInFile.open("pcf_zad4.dat");
    // avogadroInFile.open("atoms_avogadro_zad4.xyz");

    // SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
    //     beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

    // energyInFile.close();
    // posInFile.close();
    // pcfInFile.close();
    // avogadroInFile.close();

    // Zadanie piąte
    // correction = true;
    // r_start = 2.5;

    // energyInFile.open("atoms_energy_zad5.dat");
    // posInFile.open("atoms_pos_zad5.dat");
    // pcfInFile.open("pcf_zad5.dat");
    // avogadroInFile.open("atoms_avogadro_zad5.xyz");

    // SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
    //     beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

    // energyInFile.close();
    // posInFile.close();
    // pcfInFile.close();
    // avogadroInFile.close();

    // // Zadanie siódme
    // correction = true;
    // r_start = 2.5;
    // int n_start = 20;
    // int n_stop = 100;

    // std::vector<std::string> fileNames;
    // for(int n = n_start; n <= n_stop; n++){
    //     fileNames.push_back("atoms_avogadro_n_" + std::to_string(n) + ".xyz");
    //     // std::cout<<fileNames[n-n_start]<<"\n";
    // }  

    // std::ofstream energy;
    // energy.open("atoms_energy_n_zad7.dat");
    // for(int n = n_start; n <= n_stop; n++){
    //     avogadroInFile.open(fileNames[n-n_start]);
    //     V_tot = SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
    //         beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);
    //     energy << n << " "<< V_tot << "\n";
    //     avogadroInFile.close();
    // }
    // energy.close();

    // Zadanie szóste
    r_start = 2.5;
    correction = true;

    for(int i = 1; i<=5; i++){
    w_r = pow(10.0,-4);
    w_phi = pow(10.0,-2)*5;
    w_theta = pow(10.0,-2)*5;
    beta_min = 1.0;
    beta_max = 100.0;
    p = 5.0;
    std::string energyName = "atoms_energy_zad6_1_" + std::to_string(i) + ".dat";
    std::string posName = "atoms_pos_zad6_1_" + std::to_string(i) + ".dat";
    std::string pcfName = "pcf_zad6_1_" + std::to_string(i) + ".dat";
    std::string avoName = "atoms_avogadro_zad6_1_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-4);
    w_phi = pow(10.0,-2)*5;
    w_theta = pow(10.0,-2)*5;
    beta_min = 20.0;
    beta_max = 80.0;
    p = 2.5;
    energyName = "atoms_energy_zad6_2_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_2_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_2_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_2_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-4)*5;
    w_phi = pow(10.0,-2)*2;
    w_theta = pow(10.0,-2)*2;
    beta_min = 1.0;
    beta_max = 100.0;
    p = 2.0;
    energyName = "atoms_energy_zad6_3_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_3_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_3_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_3_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-4)*5;
    w_phi = pow(10.0,-2)*8;
    w_theta = pow(10.0,-2)*8;
    beta_min = 30.0;
    beta_max = 100.0;
    p = 1;
    energyName = "atoms_energy_zad6_4_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_4_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_4_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_4_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-5);
    w_phi = pow(10.0,-2)*5;
    w_theta = pow(10.0,-2)*5;
    beta_min = 1.0;
    beta_max = 10.0;
    p = 5.0;
    energyName = "atoms_energy_zad6_5_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_5_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_5_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_5_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-4)*5;
    w_phi = pow(10.0,-2)*25;
    w_theta = pow(10.0,-2)*25;
    beta_min = 1.0;
    beta_max = 100.0;
    p = 2.0;
    energyName = "atoms_energy_zad6_6_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_6_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_6_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_6_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-4)*2;
    w_phi = pow(10.0,-2)*10;
    w_theta = pow(10.0,-2)*10;
    beta_min = 25.0;
    beta_max = 75.0;
    p = 1.3;
    energyName = "atoms_energy_zad6_7_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_7_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_7_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_7_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-4)*1.8;
    w_phi = pow(10.0,-2)*5.4;
    w_theta = pow(10.0,-2)*5.4;
    beta_min = 1.0;
    beta_max = 125.0;
    p = 2.2;
    energyName = "atoms_energy_zad6_8_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_8_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_8_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_8_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-4)*1.8;
    w_phi = pow(10.0,-2)*5.4;
    w_theta = pow(10.0,-2)*5.4;
    beta_min = 1.0;
    beta_max = 200.0;
    p = 2.5;
    energyName = "atoms_energy_zad6_9_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_9_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_9_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_9_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-4);
    w_phi = pow(10.0,-2)*5;
    w_theta = pow(10.0,-2)*5;
    beta_min = 1.0;
    beta_max = 100.0;
    p = 2.0;
    energyName = "atoms_energy_zad6_10_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_10_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_10_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_10_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();

    w_r = pow(10.0,-4)*0.5;
    w_phi = pow(10.0,-2)*0.25;
    w_theta = pow(10.0,-2)*0.25;
    beta_min = 1.0;
    beta_max = 100.0;
    p = 2.0;
    energyName = "atoms_energy_zad6_11_" + std::to_string(i) + ".dat";
    posName = "atoms_pos_zad6_11_" + std::to_string(i) + ".dat";
    pcfName = "pcf_zad6_11_" + std::to_string(i) + ".dat";
    avoName = "atoms_avogadro_zad6_11_" + std::to_string(i) + ".dat";

        energyInFile.open(energyName);
        posInFile.open(posName);
        pcfInFile.open(pcfName);
        avogadroInFile.open(avoName);

        SA(n, it_max, M, w_r, w_phi, w_theta, W_all, r_start,
            beta_min, beta_max, p, correction, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

        energyInFile.close();
        posInFile.close();
        pcfInFile.close();
        avogadroInFile.close();
    }

    time(&end);
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(5); 
    std::cout << " sec \n"; 

    return 0;
}