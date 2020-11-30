//
//	Savitzsky_Golay.hpp
//  GlowCurveAnalsys
//
//  Created by Jack YuUROP 2020 Fall
//

#ifndef Savitzsky_Golay_hpp
#define Savitzsky_Golay_hpp

#include <vector>
#include <iostream>
#include <math.h>

using namespace std;

//use Savitzsky Golay algorithm to smooth data
void SG_smooth(vector<double>& y, int window_size, int size);
void transpose(vector<vector<double>> const& A, vector<vector<double>>& B, int n, int m);
vector<vector<double>> multiply(vector<vector<double>> const& A, vector<vector<double>> const& B);
double determinant(vector<vector<double>>& A, int size);
void cofactor(vector<vector<double>>& A, vector<vector<double>>& temp, int p, int q, int n);
void adjoint(vector<vector<double>>& A, vector<vector<double>>& adj);
void invert(vector<vector<double>>& A, bool neg);


#endif /* Savitzsky_Golay_hpp */