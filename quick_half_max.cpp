//
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
		auto loc = find(temp.begin(), temp.end(), TM);
		int max_index = int(loc - temp.begin());
        double energy = param[i][0];
        double dist = (c * k * TM * TM) / (energy + 2 * b * k * TM);
        double TL_index = max_index - dist;
        double TR_index = max_index + dist;
        param[i].push_back(TL_index);
        param[i].push_back(max_index);
        param[i].push_back(TR_index);
	}
}