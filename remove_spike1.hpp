//
//	remove_spike.hpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu, Jack Thieson UROP 2020 Fall
//

#ifndef remove_spike_hpp
#define remove_spike_hpp

#include <vector>
#include <algorithm>

using namespace std;

//remove the spike in raw data
void spike_elim(vector<double>& x, vector<double>& y, int span, double c, int upper, int lower);
void elim_helper(vector<double>& y, int span, double c, int index, int& process);
void elim_helper2(vector<double>& y, int span, double c, int index, int& process);

#endif /* remove_spike_hpp */