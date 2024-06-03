#include "functions.h"
#include "DSMC_2D.h"


int main(){

    srand( time( NULL ) );

    time_t start, end;
    time(&start);

    DSMC_2D ob;
    ob.read("i.dat"); // wczytujemy dane zpliku wejściowego
    ob.init(); // automatyczna inicjalizacja położeń i prędkości
    ob.write_position_velocity ("rv.dat"); // zapis ustawień początkowych
    ob.nthreads = 8; // obliczenia na jednym rdzeniu
    ob.icol = 1; // cząstki zderzają się
    ob.evolution (0.0 ,20000); // wykonujemy 20 tysięcy kroków ( tmax - nieznany )
    ob.hist_velocity_all ("hist2.dat", 5.0, 50); // zapis histogramu prędkosci do pliku
    ob.write_position_velocity ("rv.dat"); // zapis położeń i prędkości końcowych do pliku
    
    // std::string file_name = "TEST.dat";
    // file_name = "reflectance_sym_" + "TEST" + ".dat";
    // save_data(file_name, refl);

    time(&end);
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(5); 
    std::cout << " sec \n"; 

    return 0;
}