//
//  FOKModel.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 7/15/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "FOKModel.hpp"

void FOKModel(std::vector<double>& x, std::vector<double>& peak, double Tm, double Im, double E){
    //add 273.15 to Tm
    double T=0.0;
    double K = .000086173303;
    double I_t = 0.0;
    Tm += 273.15;
    double dm = (2.0*K*(Tm))/E;
    //add 273.15
    int count = 0;
    for(auto i = x.begin(); i != x.end(); ++i){
        T = *i + 273.15;
        I_t = Im*exp(1.0 +(E/(K*T))*((T-Tm)/Tm)-((T*T)/(Tm*Tm))*exp((E/(K*T))*((T-Tm)/Tm))*(1.0-((2.0*K*T)/E))-dm);
        peak[count++] = (I_t);
        
    }
}

