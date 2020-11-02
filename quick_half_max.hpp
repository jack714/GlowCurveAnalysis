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
void find_index(vector<double>& temp, vector<vector<double>>& peakParams);

#endif /* quick_half_max_hpp */