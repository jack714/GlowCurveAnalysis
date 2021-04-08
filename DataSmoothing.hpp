//
//  DataSmoothing.hpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 7/15/19.
//

#ifndef DataSmoothing_hpp
#define DataSmoothing_hpp
#include <vector>
#include <numeric>
#include <stdio.h>
double average( std::vector<double>::iterator first, std::vector<double>::iterator last, int size);
void dataSmooth(std::vector<double>& x, std::vector<double>& y);
#endif /* DataSmoothing_hpp */
