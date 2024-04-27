#include "functions_2.h"


int main(){

    srand( time( NULL ) );

    time_t start, end;
    time(&start);

    double V_tot = 0;

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
        V_tot += Brenner_i_energy(fullerene, i)/2.0;
    }

    std::cout<<"Theoretical fullerene total energy: "<<V_tot<<" eV \n";

    std::ofstream energyInFile, posInFile, pcfInFile, avogadroInFile;
    energyInFile.open("atoms_energy_01.dat");
    posInFile.open("atoms_pos_01.dat");
    pcfInFile.open("pcf_01.dat");
    avogadroInFile.open("atoms_avogadro_01.xyz");

    SA(60, pow(10,5), 100, pow(10,-4), pow(10,-2)*5, pow(10,-2)*5, pow(10,-4), 3.5, 1.0, 100, 2, &energyInFile, &posInFile, &pcfInFile, &avogadroInFile);

    energyInFile.close();
    posInFile.close();
    pcfInFile.close();
    avogadroInFile.close();

    time(&end);
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(5); 
    std::cout << " sec \n"; 

    return 0;
}