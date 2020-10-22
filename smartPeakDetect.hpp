//
//  smartPeakDetect.hpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 7/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#ifndef smartPeakDetect_hpp
#define smartPeakDetect_hpp

#include <string>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <locale>
#include <numeric>
#include <fstream>
#include <sstream>
#include "FOKModel.hpp"
#include "DataSmoothing.hpp"

void smartPoints(std::vector<double>& x, std::vector<double>& y, std::vector<int>& minimum,std::vector<int>& maxima,std::vector<double> derivative,std::vector<double> secDerivative,std::vector<int>& inflectPnt);

void pointsParams(std::vector<double>& x, std::vector<double>& y, std::vector<int>&maxima, std::vector<int>& minima, std::vector<std::vector<double>>& peakParams);


//void dataSmooth(std::vector<double>& x, std::vector<double>& y, std::vector<double>& xNew, std::vector<double>& yNew);
double activation(double TL, double TR, double TM);

// Sudo Main function for smartPeakDetect, take in temperature and count data then record peaks in peakParams
void findPeaks(std::vector<double>& x, std::vector<double>& y, std::vector<std::vector<double>>& peakPrams, std::string output_dir);

//double average(std::vector<double>::iterator first, std::vector<double>::iterator last, int size);

void firstDeriv(std::vector<double>& x, std::vector<double>& y, std::vector<double>& derivative);

void nonMaxPeaks(std::vector<double>& x, std::vector<double>& y, std::vector<double> secDerivative, std::vector<int>& maxima, std::vector<int>& minima, std::vector<std::vector<double>>& peakParams, std::string output_dir);

void secDeriv(std::vector<double>& x, std::vector<double>& y, std::vector<double>& derivative);

void printFindings(std::vector<double>& x, std::vector<double>& y, std::vector<int>& minimum,std::vector<int>& maxima, std::vector<int>& inflectPnt, std::string dir);

void write(std::vector<std::vector<double>> glow_curves,std::vector<double> y,std::vector<double> x, std::string output_name);
#endif /* smartPeakDetect_hpp */
