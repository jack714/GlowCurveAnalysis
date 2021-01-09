//
//	background_remove.cpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu, Jack T, UROP 2020 Fall
//
#include "background_remove.hpp"
#include "smartPeakDetect.hpp"


void remove_back(vector<double>& x, vector<double>& y) {
    vector<double> firstDir(x.size(), 0.0);
    //calculate first derivative of count data using function from SmartPeakDetect
    firstDeriv(x, y, firstDir);
    //derivative change at the left and right to idenify the start of actual data
    double left = 0.02;
    double right = -.01;
    int leftPoint = 0;
    int rightPoint = 0;
    //find the index of the point where the derivative satisfies left cap
    for (int i = 0; i < static_cast<int>(firstDir.size()); i++) {
        if (firstDir[i] >= left) {
            leftPoint = i;
            break;
        }
    }
    //find the index of the point where the derivative satisfies right cap
    for (int i = static_cast<int>(firstDir.size()) - 1; i >= 0; i--) {
        if (firstDir[i] >= right) {
            rightPoint = i;
            break;
        }
    }
    //divide the start of data to the leftPoint in half and take each's count average
    int middle = leftPoint / 2;
    double firstCount = 0.0;
    double secCount = 0.0;
    for (int i = 0; i < middle; i++)
        firstCount += y[i];
    firstCount /= middle;
    for (int i = middle; i < leftPoint; i++)
        secCount += y[i];
    secCount /= middle;
    //take the mid temperature over the two half range
    double firstTemp = x[leftPoint / 4];
    double secTemp = x[leftPoint / 4 * 3];
    //cout << firstTemp << " " << secTemp << endl;
    //cout << firstCount << " " << secCount << endl;
    //approximate a line from the two
    double slope = (secCount - firstCount) / (secTemp - firstTemp);
    double c = firstCount - slope * firstTemp;
    //check if the right part approximation is strictly lower than the observed data
    bool check = true;
    for (int i = rightPoint; i < static_cast<int>(x.size()); i++) {
        if (slope * x[i] + c > y[i]) {
            check = false;
            break;
        }
    }
    int size = static_cast<int>(x.size());
    if (!check) {
        double leftCount = (firstCount + secCount) / 2;
        double leftTemp = x[middle];
        double rightCount = 0.0;
        double rightTemp = 0.0;
        for (int i = rightPoint; i < size; i++)
            rightCount += y[i];
        rightCount /= (size - rightPoint);
        rightTemp = x[(size - rightPoint) / 2];
        slope = (rightCount - leftCount) / (rightTemp - leftTemp);
        c = leftCount - slope * leftTemp;
    }
    for (int i = 0; i < size; i++) 
        y[i] -= slope * x[i] + c;
    vector<double> xTemp;
    vector<double> yTemp;
    for (int i = leftPoint; i < rightPoint; i++) {
        xTemp.push_back(x[i]);
        yTemp.push_back(y[i]);
    }
    swap(x, xTemp);
    swap(y, yTemp);
    //cout << slope << " " << c;
}