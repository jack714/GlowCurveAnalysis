//
//	remove_spike.cpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu, Jack Thieson UROP 2020 Fall
//

#include "remove_spike.hpp"

void spike_elim(vector<double>& x, vector<double>& y, int span, double c) {
    vector<int> alter;
    //process all data points and store the left point of the left boundary point to the vector alter
    for (int i = 0; i < static_cast<int>(x.size()) - span; i++) {
        int process = -1;
        elim_helper1(y, span, c, i, process);
    }
    for (int i = 0; i < static_cast<int>(x.size()) - span; i++) {
        int process = -1;
        elim_helper2(y, span, c, i, process);
        if (process != -1 && process + span < static_cast<int>(y.size()) - 1) {
            alter.push_back(process);
            if (process + span < static_cast<int>(y.size()) - span) {
                for (int j = 1; j < span; j++) {
                    alter.push_back(process + j);
                }
            }
        }
    }
    //if alter is not empty then continue to further process
    int iteration = 0;
    while (!alter.empty() && iteration < 4) {
        //change the parameter
        span += 2;
        c -= 0.01;
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
            elim_helper2(y, span, c, i, process);
        }
        iteration++;
    }
    for (int i = 0; i < static_cast<int>(x.size()) - span; i++) {
        int process = -1;
        elim_helper2(y, 3, c, i, process);
    }
    for (int i = 0; i < static_cast<int>(x.size()) - span; i++) {
        int process = -1;
        elim_helper2(y, 5, c, i, process);
    }
}

void elim_helper1(vector<double>& y, int span, double c, int index, int& process) {
    if ((y[index] + y[index + span - 1]) / 2 != 0 && y[index + span / 2] > 3) {
        //double d = 0.0;
        //for (int i = index; i < index + span; i++)
        //    d += y[i];
        //double avg = d / span;
        if ((y[index + 1] < 1 / c * (y[index] + y[index + 2]) / 2) && (y[index + 2] > c* (y[index + 1] + y[index + 3]) / 2)) {
            double change = (y[index + 3] - y[index]) / 3;
            y[index + 1] = y[index] + change;
            y[index + 2] = y[index] + 2 * change;
            process = index - 1;
        }
    }
}
void elim_helper2(vector<double>& y, int span, double c, int index, int& process) {
    if ((y[index] + y[index + span - 1]) / 2 != 0 && y[index + span / 2] > 3) {
        //double d = 0.0;
        //for (int i = index; i < index + span; i++)
        //    d += y[i];
        //double avg = d / span;
        //if ((y[index + 1] < 1 / c * (y[index] + y[index + 2]) / 2) && (y[index + 2] > c* (y[index + 1] + y[index + 3]) / 2)) {
        //    double change = (y[index + 3] - y[index]) / 3;
        //    y[index + 1] = y[index] + change;
        //    y[index + 2] = y[index] + 2 * change;
        //    process = index - 1;
        //}
        //if ((y[index + span / 2] > c* (y[index] + y[index + span - 1]) / 2) || (y[index + span / 2] < 1 / c * (y[index] + y[index + span - 1]) / 2)) {
        if ((y[index + span / 2] > c* (y[index] + y[index + span - 1]) / 2)) {
            //if((y[index + span / 2] > c * avg)) {
                //double change = (y[index + span - 1] - y[index]) / (span - 1);
                //for (int j = index + 1; j < index + span - 1; j++) {
                //    y[j] = y[j - 1] + change;
                //}
                //y[index + span / 2] =  c * (y[index] + y[index + span - 1]) / 2;
            y[index + span / 2] = (y[index] + y[index + span - 1]) / 2;
            //y[index + span / 2] = avg;
            process = index - 1;
        }
        else if (y[index + span / 2] < 1 / c * (y[index] + y[index + span - 1]) / 2) {
        //else if (y[index + span / 2] < 1 / c * avg) {
            y[index + span / 2] = 1 / c * (y[index] + y[index + span - 1]) / 2;
        //    y[index + span / 2] = (y[index] + y[index + span - 1]) / 2;
            //y[index + span / 2] = avg;
        //    process = index - 1;
        }
    }
}
