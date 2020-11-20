//
//	remove_spike.cpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu, Jack Thieson UROP 2020 Fall
//

#include "remove_spike.hpp"
#include "smartPeakDetect.hpp"

void spike_elim(vector<double>& x, vector<double>& y) {
    //make sure only process minimas that are separated within dist, so that not treating actual peak as spike
    int dist = 10;
    //constant used to qualify if a point is a spike
    double c = 2;
	vector<double> derivative(x.size(), 0.0);
	//call firstDeriv function from smartPeakDetect.cpp and populate derivative vector with first derivative data
	firstDeriv(x, y, derivative);
	int size = int(derivative.size());
	vector<int> minimum;
    //find all the local minimums
    for (int i = 1; i < size; i++)
    {
        //minimum is recorded as first derivative change from negative to positive
        if (derivative[i] > 0.0 && derivative[i - 1] < 0.0)
        {
            minimum.push_back(i);
        }
    }
    minimum.push_back(int(x.size()) - 1);
    for (int j = 0; j < static_cast<int>(minimum.size()) - 1; j++) {
        //only process data that are within dist apart
        if (minimum[j + 1] - minimum[j] < dist) {
            int mid = minimum[j] + (minimum[j + 1] - minimum[j]) / 2;
            //recognize spike if the y count is larger than c times the average of the two minima points
            if (y[mid] > c* (y[j] + y[j + 1]) / 2) {
                //use straight line approximation to remove spike
                int span = minimum[j + 1] - minimum[j];
                double change = (y[j + 1] - y[j]) / span;
                for (int k = minimum[j] + 1; k < minimum[j + 1]; k++) {
                    y[k] = y[k - 1] + change;
                }
            }
        }
    }
}

void spike_elim2(vector<double>& x, vector<double>& y) {
    double c = 2;
    for (int i = 0; i < static_cast<int>(x.size()) - 10; i++) {
        if (y[i + 5] > c* (y[i] + y[i + 9]) / 2) {
            double change = (y[i + 9] - y[i]) / 9);
            for (int j = i + 1; j < i + 9; j++) {
                y[j] = y[j - 1] + change;
            }
        }
    }
}