//
//  OTORModel.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 7/15/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include <stdio.h>
#include "OTORModel.hpp"


void OTORModel(std::vector<double>& x, std::vector<double>& peak, double Tm, double Im, double E){
    //add 273.15 to Tm
    double T=0.0;
    const double K = .000086173303;
    const double c = 2.99792458e8;
    const double R = 8.314462618;
    double I_t = 0.0;
    Tm += 273.15;
    int count = 0;
    double integral = 0.0;
    double d = -E/(K*Tm);
    double z1m = (1/c) - log(c) + (E*exp(-d)) / ((K*Tm*Tm)*(-1.05 * (pow(R,(1.26))) + 1));
    z1m *= Tm * exp(d) + (E/K) * expint(d);
    for(auto i = x.begin(); i != x.end(); ++i){
        T = *i + 273.15;
        d = -E/(K*T);
        integral += T * exp(d) + (E/K) * expint(d);
        double z1 = (1/c) - log(c) + (E*exp(-d)) / ((K*Tm*Tm)*(-1.05 * pow(R,(1.26)) + 1));
        z1 *= integral;
        
        I_t = Im * (lambert_w0(exp(z1m)) + pow(lambert_w0(exp(z1m)),2.0))/(lambert_w0(exp(z1)) + pow(lambert_w0(exp(z1)),2.0));
        double temp = (-E/K) * ((1/T)-(1/Tm));
        I_t *= exp(temp);
        peak[count++] = I_t;
    }
}
