#include "utils.h"


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
