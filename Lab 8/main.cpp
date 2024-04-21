#include "functions.h"



int main(){

    std::fstream pos;
    pos.open("pos.dat", std::fstream::out);


    pos.close();

    return 0;

}