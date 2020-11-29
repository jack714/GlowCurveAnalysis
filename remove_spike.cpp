//
//	remove_spike.cpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu, Jack Thieson UROP 2020 Fall
//

#include "remove_spike.hpp"

void spike_elim(vector<double> & x, vector<double> & y, int span, double c) {
    vector<int> alter;
    //process all data points and store the left point of the left boundary point to the vector alter
    for (int i = 0; i < static_cast<int>(x.size()) - span; i++) {
        int process = -1;
        elim_helper(y, span, c, i, process);
        if (process != -1 && process + span < static_cast<int>(y.size()) - 1) {
            alter.push_back(process);
            if (process + span < static_cast<int>(y.size()) - span)
                for (int j = 1; j < span; j++) {
                    alter.push_back(process + j);
                }
        }
    }
    //if alter is not empty then continue to further process
    int iteration = 0;
    while (!alter.empty() && iteration < 4) {
        //change the parameter
        span += 2;
        c -= 0.05;
        //remove points that will go out of range
        int j = 0;
        while (j < static_cast<int>(alter.size())) {
            if (alter[j] + span > static_cast<int>(y.size()))
                alter.erase(alter.begin() + j);
            else
                j++;
        }
        //only process the boundary points stored in alter
        for (int i : alter) {
            int process = -1;
            elim_helper(y, span, c, i, process);
        }
        iteration++;
    }
}

void elim_helper(vector<double>& y, int span, double c, int index, int& process) {
    if ((y[index] + y[index + span - 1]) / 2 != 0 && y[index + span / 2] > 3) {
        if ((y[index + span / 2] > c* (y[index] + y[index + span - 1]) / 2) || (y[index + span / 2] < 1 / c * (y[index] + y[index + span - 1]) / 2)) {
            double change = (y[index + span - 1] - y[index]) / (span - 1);
            for (int j = index + 1; j < index + span - 1; j++) {
                y[j] = y[j - 1] + change;
            }
            process = index - 1;
        }
    }
}
