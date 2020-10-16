//
//  DataSmoothing.hpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 7/15/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#ifndef DataSmoothing_hpp
#define DataSmoothing_hpp
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdio.h>
double average( std::vector<double>::iterator first, std::vector<double>::iterator last, int size);
void dataSmooth(std::vector<double>& x, std::vector<double>& y);
#endif /* DataSmoothing_hpp */
