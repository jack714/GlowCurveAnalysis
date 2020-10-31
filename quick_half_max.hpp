//
//  GlowCurveAnalsys
//
//  Created by Jack Yu UROP 2020 Fall
//

#ifndef quick_half_max_hpp
#define quick_half_max_hpp

#include <stdio.h>
#include <string>
#include <cmath>
#include <vector>
#include <locale>
#include <numeric>
#include <fstream>
#include <algorithm>
#include "FOKModel.hpp"

using namespace std;

//find the half_max points for every peak and populate peakParams accordingly
void find_index(vector<double>& temp, vector<double>& count, vector<vector<double>>& peakParams);

// activation formula used in find_half_max helper function
double energy(double TL, double TR, double TM);

#endif /* quick_half_max_hpp */