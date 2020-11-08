//
//  background_subtraction.cpp
//  GlowCurveAnalsys
//
//  Created by jack thiesen on 11/4/2020.
//

#include "background_subtraction.hpp"

void bg_subtract(std::vector<double>& x, std::vector<double>& y) {
	// while under 90 degrees C
	int i = 0;
	double avg_start = 0;
	while (x[i] < 90) {
		// and the counts are low (preventing spikes) (USE AVERAGE LATER!!!)
		if (y[i] < 5)
			avg_start += y[i];
		++i;
	}
	// average the counts
	if(i != 0)
		avg_start /= i;
	
	// perform subtraction
	for(size_t j = 0; j < y.size(); j++) {
		y[j] -= avg_start;
		// if negative set to 0
		if (y[j] < 0)
			y[j] = 0.0;
	}
	/*
	Below has problems at low dose, fix later using averages probably
	
	// iterate backwards averaging counts until the abs(slope) is greater than 2
	double avg_end = x[x.length() - 1];
	i = x.length() - 1;
	while ( y[i] < 15 && abs((y[i] - y[i - 1])/2) < 2) {
		
		--i;
	}
	// if avg_end fails to get a value or has too high of a value, assume it to be 6
	if (avg_end > 10.0)
		avg_end = 6.0;
	// create slope of line to be subtracted
	
	// perform the subtraction, set any negative numbers to 0
	for(j = 0; j < y.length(); j++) {
		y[j] -= ;
		if (y[j] < 0)
			y[j] = 0.0;
	}
	*/
}