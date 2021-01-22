//
//	remove_spike.cpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu, Jack Thieson UROP 2020 Fall
//

#include "remove_spike.hpp"

void spike_elim(vector<double>& x, vector<double>& y, int span, double c) {
    //vector<int> alter;
    //process all data points and store the left point of the left boundary point to the vector alter
    for (int i = 0; i < static_cast<int>(x.size()) - span; i++) {
        //int process = -1;
        elim_helper(y, span, c, i);
        //if (process != -1 && process + span < static_cast<int>(y.size()) - 1) {
        //    alter.push_back(process);
        //    if(process + span < static_cast<int>(y.size()) - span)
        //    for (int j = 1; j < span; j++) {
        //        alter.push_back(process + j);
        //    }
        //}
    }
    //if alter is not empty then continue to further process
    int iteration = 0;
    //while (!alter.empty() && iteration < lower) {
    //while(iteration < lower) {
        //change the parameter
    //    if(span < 15)
    //        span += 2;
    //    c -= 0.05;
        //remove points that will go out of range
        //int j = 0;
        //while (j < static_cast<int>(alter.size())) {
        //    if (alter[j] + span > static_cast<int>(y.size()))
        //        alter.erase(alter.begin() + j);
        //    else
        //        j++;
        //}
        //only process the boundary points stored in alter
        //for (int i : alter) {
        //    int process = -1;
    //    int process = -1;
    //    for(int i = 0; i < static_cast<int>(y.size()) - span; i++)
    //        elim_helper(y, span, c, i, process);
        //}

    //    iteration++;
    //}
    //iteration = 0;
    //while (!alter.empty() && iteration < upper) {
    while(iteration < 5) {
        //change the parameter
        span += 2;
        c -= 0.025;
        //remove points that will go out of range
        //int j = 0;
        //while (j < static_cast<int>(alter.size())) {
        //    if (alter[j] + span > static_cast<int>(y.size()))
        //        alter.erase(alter.begin() + j);
        //    else
        //        j++;
        //}
        //only process the boundary points stored in alter
        //for (int i : alter) {
        //    int process = -1;
        //    elim_helper2(y, span, c, i, process);
        //}
        for (int i = 0; i < int(y.size()) - span; i++)
            elim_helper(y, span, c, i);
        iteration++;
    }
}

void elim_helper(vector<double>& y, int span, double c, int index) {
    if ((y[index] + y[index + span - 1]) / 2 != 0 && y[index + span / 2] > 3) {
        //if ((y[index + span / 2] > c * (y[index] + y[index + span - 1]) / 2) || (y[index + span / 2] < 1/c * (y[index] + y[index + span - 1]) / 2)) {
        //if ((y[index + span / 2] > c * pow(y[index] * y[index + span - 1], 0.5)) || (y[index + span / 2] < 1 / c * pow(y[index] * y[index + span - 1], 0.5))) {
        if ((y[index + span / 2] > c* (y[index] + y[index + span - 1]) / 2)) {
            double change = (y[index + span - 1] - y[index]) / (span - 1);
            for (int j = index + 1; j < index + span - 1; j++) {
                y[j] = y[j - 1] + change;
            }
        }
    }
}

void elim_helper2(vector<double>& y, int span, double c, int index, int& process) {
    if ((y[index] + y[index + span - 1]) / 2 != 0 && y[index + span / 2] > 3) {
        //vector<double> data;
        //for (int i = 0; i < span; i++) {
        //    data.push_back(y[index + i]);
        //}
        //sort(data.begin(), data.end());
        //double med = data[span / 2];
        //if (y[index + span / 2] > c* med) {
        if ((y[index + span / 2] > c* (y[index] + y[index + span - 1]) / 2)) {
            double change = (y[index + span - 1] - y[index]) / (span - 1);
            for (int j = index + 1; j < index + span - 1; j++) {
                y[j] = y[j - 1] + change;
            }
            process = index - 1;
        }
    }
}

