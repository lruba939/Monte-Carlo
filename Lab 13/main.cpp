#include "functions.h"
#include "DSMC_2D.h"


int main(){

    srand( time( NULL ) );

    time_t start, end;
    time(&start);

    // DSMC_2D ob;
    // ob.read("/home/ruba/Monte-Carlo/Lab 13/Zad1/i1.dat"); // wczytujemy dane zpliku wejściowego
    // ob.init(); // automatyczna inicjalizacja położeń i prędkości
    // ob.write_position_velocity ("rv_left.dat"); // zapis ustawień początkowych
    
    // DSMC_2D ob2;
    // ob2.read("/home/ruba/Monte-Carlo/Lab 13/Zad2/i2.dat"); // wczytujemy dane zpliku wejściowego
    // ob2.init(); // automatyczna inicjalizacja położeń i prędkości
    // ob2.write_position_velocity ("rv_right.dat"); // zapis ustawień początkowych

    DSMC_2D ob3;
    ob3.read("/home/ruba/Monte-Carlo/Lab 13/i.dat");
    ob3.init();
    ob3.nthreads=4;
    ob3.icol=1;
    ob3.evolution(0.0,2000);
    ob3.write_position_velocity("rv_end.dat");


    time(&end);
    double time_taken = double(end - start); 
    std::cout << "Time taken by program is : " << std::fixed << time_taken << std::setprecision(5); 
    std::cout << " sec \n"; 

    return 0;
}