#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <bits/stdc++.h>
#include <iomanip>
#include <time.h>


double uni_rand();

void save_data(std::string outFileName, std::vector<std::vector<double>> data);

std::vector<std::vector<double>> load_data(std::string inFileName);


#endif