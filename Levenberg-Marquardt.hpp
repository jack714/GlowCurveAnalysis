//
//  First_Order_kinetics.hpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#ifndef First_Order_kinetics_hpp
#define First_Order_kinetics_hpp

#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <stdio.h>

class First_Order_Kinetics{
private:
    std::vector<double> count_data, curve, curve_areas, temp_data;
    std::vector<std::vector<double>> glow_curves, peakParams;
    const int MAX_ITER = 1000;
    double k = .000086173303;
    double totalArea;
public:
    First_Order_Kinetics(std::pair<std::vector<double>,std::vector<double>>, std::vector<std::vector<double>>);
    //double activation_energy(int TL_index,int TM_index,int TR_index);
    double glow_curve();
    std::vector<double> initial_guess(std::vector<double> &curve,int i);
//    double Func(const double input, const std::vector<double> params);
    double Func2(const double input, const std::vector<double> params);

    std::vector<std::vector<double>> jacobian(int max_index,int TL_index,int TR_index,double E,double Tm,double Im);
    void transpose(std::vector<std::vector<double>> const &A,std::vector<std::vector<double>> &B, int n, int m);
    std::vector<std::vector<double>> multiply(std::vector<std::vector<double>> const &A,std::vector<std::vector<double>> const &B);
    void invert(std::vector<std::vector<double>> &A, bool neg);
    double determinant(std::vector<std::vector<double>> &A, int size);
    void cofactor(std::vector<std::vector<double>> &A, std::vector<std::vector<double>> &temp, int p, int q, int n);
//    double Deriv(const double input, const std::vector<double> params, int n);
    double Deriv2(const double input, const std::vector<double> params, int n);

    double dotProduct(std::vector<double> A, std::vector<double> B);
    std::vector<std::vector<double>> Identity(int num, double lambda);
    std::vector<double> vec_matrix_multi(std::vector<std::vector<double>> const &A,std::vector<double> const &B);
    void LevenbergMarquardt2(const std::vector<double> &outputs, std::vector<double> &params);
    void LevenbergMarquardt(const std::vector<double> &outputs, std::vector<std::vector<double>> &params, double &FOM);
    void adjoint(std::vector<std::vector<double>> &A,std::vector<std::vector<double>> &adj);
    std::vector<std::vector<double>> jacobian(const std::vector<std::vector<double>> &inputs, std::vector<double> params);
    std::vector<std::vector<double>> return_glow_curve(){
        return glow_curves;
    }
    std::vector<double> return_curve_areas(){
        return curve_areas;
    }
    double area() {
        return totalArea;
    }
};

#endif
