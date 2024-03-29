//
//	background_remove.cpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu, Jack T, UROP 2020 Fall
//
#include "background_remove.hpp"
#include "smartPeakDetect.hpp"


vector<double> remove_back(vector<double>& x, vector<double>& y) {
    vector<double> firstDir(x.size(), 0.0);
    //vector<double> secDir(x.size(), 0.0);
    //calculate first derivative of count data using function from SmartPeakDetect
    firstDeriv(x, y, firstDir);
    //secDeriv(x, y, secDir);
    //derivative change at the left and right to idenify the start of actual data
    //old gca double left = 0.02;
    double left = 0.02;
    //double left = -.003;
    //old gca double right = -0.009;
    double right = -0.009;
    //double right = 0.001;
    int leftPoint = 2;
    int rightPoint = int(x.size()) - 2;
    //find the index of the point where the derivative satisfies left cap
    //for (int i = 0; i < static_cast<int>(firstDir.size()); i++) {
    //    if (firstDir[i] >= left) {
    //        leftPoint = i;
    //        break;
    //    }
    //}
    //find the index of the point where the derivative satisfies right cap
    //for (int i = static_cast<int>(firstDir.size()) - 1; i >= 0; i--) {
    //    if (firstDir[i] <= right) {
    //        rightPoint = i;
    //        break;
    //    }
    //}
    //divide the start of data to the leftPoint in half and take each's count average
    int middle = leftPoint / 2;
    double firstCount = 0.0;
    double secCount = 0.0;
    for (int i = 0; i < leftPoint; i++)
        firstCount += y[i];
    firstCount /= (leftPoint+1);
    for (int i = rightPoint; i < static_cast<int>(y.size()); i++)
        secCount += y[i];
    secCount /= (static_cast<int>(y.size())-rightPoint+1);
    //take the mid temperature over the two half range
    double firstTemp = x[leftPoint/2];
    double secTemp = x[(static_cast<int>(y.size()) + rightPoint + 1) / 2];
    //cout << firstTemp << " " << secTemp << endl;
    //cout << firstCount << " " << secCount << endl;
    //approximate a line from the two
    double slope = (secCount - firstCount) / (secTemp - firstTemp);
    double c = firstCount - slope * firstTemp;
    //check if the right part approximation is strictly lower than the observed data
    bool check = true;
    for (int i = rightPoint; i < static_cast<int>(x.size()); i++) {
        if (slope * x[i] + c > y[i]) {
            //cout << "new: " << slope * x[i] + c << endl;
            //cout << "x is: " << i << endl;
            //cout << "y: " << y[i] << endl;
            check = false;
            break;
        }
    }
    //cout << slope << " " << c << endl;
    int size = static_cast<int>(x.size());
    //if (!check) {
    //    //cout << "yes!" << endl;
    //    double leftCount = (firstCount + secCount) / 2;
    //    double leftTemp = x[middle];
    //    double rightCount = 0.0;
    //    double rightTemp = 0.0;
    //    for (int i = rightPoint; i < size; i++)
    //        rightCount += y[i];
    //    rightCount /= (size - rightPoint);
    //    rightTemp = x[rightPoint + (size - rightPoint) / 2];
    //    slope = (rightCount - leftCount) / (rightTemp - leftTemp);
    //    c = leftCount - slope * leftTemp;
    //    //cout << leftTemp << " " << leftCount << endl;
    //    //cout << rightTemp << " " << rightCount << endl;
    //}
    for (int i = 0; i < size; i++) {
        y[i] -= slope * x[i] + c;
        if (y[i] < 0)
            y[i] = 0;
    }
    //vector<double> xTemp;
    //vector<double> yTemp;
    //for (int i = leftPoint; i < rightPoint; i++) {
    //    xTemp.push_back(x[i]);
    //    yTemp.push_back(y[i]);
    //}
    //swap(x, xTemp);
    //swap(y, yTemp);
    //cout << slope << " " << c;
    //cout << xTemp[rightPoint];
    return firstDir;
}
