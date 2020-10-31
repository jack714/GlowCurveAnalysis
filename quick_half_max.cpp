//
//  GlowCurveAnalsys
//
//  Created by Jack Yu UROP 2020 Fall
//

#include "quick_half_max.hpp"

//find the half_max points for every peak and populate peakParams accordingly
void find_index(vector<double>& temp, vector<double>& count, vector<vector<double>>& param) {
	//for each find its index and then find half_max index
	for (int i = 0; i < static_cast<int>(param.size()); i++) {
		//find the index of the peak
		auto loc = find(temp.begin(), temp.end(), param[i][1]);
		int max_index = int(loc - temp.begin());
		double height = count[max_index];
		//find the right half_max peak's index
		auto right_half_loc = find(temp.begin() + max_index, temp.end(), 0.5 * height);
		int right_half = right_half_loc - temp.begin();
		//using symmetry to find the left half_max index;
		int left_half = max_index - (right_half - max_index);
		//push all three indexes to peakParams
        param[i].push_back(left_half);
        param[i].push_back(max_index);
        param[i].push_back(right_half);
		//calculate activation energy for the peaks
        param[i][0] = energy(temp[left_half], temp[max_index], temp[right_half]);
	}
}

// activation formula used in find_half_max helper function
double energy(double TL, double TR, double TM)
{
    double m_g = 0.0, E = 0.0, C = 0.0, K = .000086173303;
    double b = 0.0, t = 0.0, d = 0.0, w = 0.0;
    TL += 273.15;
    TM += 273.15;
    TR += 273.15;

    t = TM - TL;
    if (t == 0)
    {
        t = 1;
    }

    w = TR - TL;
    d = TR - TM;
    m_g = d / w;
    b = 1.58 + 4.2 * (m_g - 0.42);
    C = 1.51 + (3 * (m_g - 0.42));

    E = ((C * K * (TM * TM)) / t) - (b * (2 * K * TM));
    if (E > 3.0)
    {
        E = 3;
    }
    return E;
}