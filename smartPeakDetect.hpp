//
//  smartPeakDetect.hpp
//  GlowCurveAnalsys
//
//  Initially created by jeremy hepker on 7/9/19.
//
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//

#ifndef smartPeakDetect_hpp
#define smartPeakDetect_hpp

#include <vector>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "FOKModel.hpp"
#include "DataSmoothing.hpp"
using namespace std;

void smartPoints(std::vector<double>& x, std::vector<double>& y, std::vector<int>& minimum,std::vector<int>& maxima,std::vector<double> derivative,std::vector<double> secDerivative,std::vector<int>& inflectPnt);

//half width half max method, populate peakParams with activation data,
// and temperature, count of half width half max, along with its left, middle, and right index
void pointsParams(std::vector<double>& x, std::vector<double>& y, std::vector<int>&maxima, std::vector<int>& minima, std::vector<std::vector<double>>& peakParams);
    
//void dataSmooth(std::vector<double>& x, std::vector<double>& y, std::vector<double>& xNew, std::vector<double>& yNew);
double activation(double TL, double TR, double TM);

//helper function that runs half width half max method once and returns TM_index
int find_half_max(int index, std::vector<double>& x,std::vector<double>& y, std::vector<int>& maxima, std::vector<int>& minima, std::vector<std::vector<double>>& peakParams);

// Sudo Main function for smartPeakDetect, take in temperature and count data then record peaks in peakParams
std::vector<double> findPeaks(std::vector<double>& x, std::vector<double>& y, std::vector<std::vector<double>>& peakPrams, std::string output_dir, std::vector<double>& firstDir);

//double average(std::vector<double>::iterator first, std::vector<double>::iterator last, int size);

void firstDeriv(std::vector<double>& x, std::vector<double>& y, std::vector<double>& derivative);

void nonMaxPeaks(std::vector<double>& x, std::vector<double>& y, std::vector<double> secDerivative, std::vector<int>& maxima, std::vector<int>& minima, std::vector<std::vector<double>>& peakParams, std::string output_dir);

void secDeriv(std::vector<double>& x, std::vector<double>& y, std::vector<double>& derivative);

void printFindings(std::vector<double>& x, std::vector<double>& y, std::vector<int>& minimum,std::vector<int>& maxima, std::vector<int>& inflectPnt, std::string dir);

void write(std::vector<std::vector<double>> glow_curves,std::vector<double> y,std::vector<double> x, std::string output_name);
#endif /* smartPeakDetect_hpp */
