//
//  First_Order_kinetics.hpp
//  GlowCurveAnalsys
//
//  Initially created by jeremy hepker on 1/27/19.
//
//  Modified and re-organized by Jack Yu UROP 2020 Fall

#ifndef First_Order_kinetics_hpp
#define First_Order_kinetics_hpp

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>

class First_Order_Kinetics{
private:
    std::vector<double> count_data, curve, curve_areas, temp_data, orig_sig_deriv, decon_sig_deriv;
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
    //calculate first derivative
    void deriv(std::vector<double>& x, std::vector<double>& y, std::vector<double>& derivative);

    //populate decon_sig_deriv for the sum of the deconvolute curve
    void update_deriv();

    //use orig_sig_deriv and decon_sig_deriv
    //(osd[i] - dsd[i]) / (sum of osd) * 100 ----> FOM in %
    void deriv_FOM();
};

#endif
