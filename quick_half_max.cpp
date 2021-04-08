//
//  quick_half_max.cpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu UROP 2020 Fall
//

#include "quick_half_max.hpp"

//find the half_max points for every peak and populate peakParams accordingly
void find_index(vector<double>& temp, vector<vector<double>>& param) {
	//for each peak derive the left anbd right half_max points from activation energy
    //assume the distance from t_max to t_left and t_right is the same
    //formula: E = ((C * K * (TM * TM)) / t) - (b * (2 * K * TM));
    double b = 1.916;
    double c = 1.75;
    double k = 0.000086173303;
	for (int i = 0; i < static_cast<int>(param.size()); i++) {
		//find the index of the peak
        double TM = param[i][1];
        int max = 0;
        while (temp[max] - TM < 0.1 && max < static_cast<int>(temp.size()) - 1) {
            max++;
        }
        int max_index = max;
		//auto loc = find(temp.begin(), temp.end(), TM);
        TM += 273.15;
		//int max_index = int(loc - temp.begin());
        double energy = param[i][0];
        double dist = (c * k * TM * TM) / (energy + 2 * b * k * TM);
        TM -= 273.15;
        double left_temp = TM - dist;
        int index = max_index;
        while (temp[index] - left_temp > 0.1 && index > 1) {
            index--;
        }
        int left_index = index;
        int right_index = abs(left_index - max_index) + max_index;
        if (right_index >= static_cast<int>(temp.size()) - 1)
            right_index = static_cast<int>(temp.size()) - 3;
        param[i][3] = left_index;
        param[i][4] = max_index;
        param[i][5] = right_index;
	}
}