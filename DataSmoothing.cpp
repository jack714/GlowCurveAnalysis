//
//  DataSmoothing.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 7/15/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "DataSmoothing.hpp"
void dataSmooth(std::vector<double>& x, std::vector<double>& y){
    //Smooth data by averaging plateau'ed temperatures
    int first = 0, sec = 0;
    std::vector<double> yNew, xNew;
    xNew.reserve(x.size());
    yNew.reserve(y.size());
    for(int i = 0; i < int(x.size());++i){
        if(x[i] == x[first]){
            ++sec;
        }else{
            xNew.push_back(x[i]);
            yNew.push_back(average(y.begin()+first, y.begin()+sec, (sec - first)));
            first = i;
            sec = i+1;
        }
    }
    for(int j = 0; j < 2; ++j){
        for(int i = 2; i < int(yNew.size()-2); ++i){
            yNew[i] = (yNew[i-2] + 2.0 * yNew[i-1] + 3.0 * yNew[i] + 2.0 * yNew[i+1] + yNew[i+2])/9;
        }
    }
    x = xNew;
    y = yNew;
}
double average( std::vector<double>::iterator first, std::vector<double>::iterator last, int size){
    double sum = accumulate(first, last, 0.0);
    sum /= double(size);
    return sum;
}
