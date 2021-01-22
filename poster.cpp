//
//  main.cpp
//  GlowCurveAnalsys
//
//  Initially created by jeremy hepker on 1/9/19.
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//

#include <iostream>
#include <getopt.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h> 
#include <random>
#include <math.h>
#include "remove_spike.hpp"
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
//if unistd.h can't be included then comment out the else clause
//#else
//#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;


double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

int main() {
    ofstream file1;
    string path = "C:/Users/jack0/Desktop/output.csv";
    file1.open(path);
    double lower_bound = 0;
    double upper_bound = 10;
    //std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    //std::default_random_engine re;
    //double e = unif(re);
    //double g = unif(re);
    double pi = 2 * acos(0.0);
    vector<double> x;
    vector<double> y;
    bool mode = true;
    if (mode) {
        for (double i = 0; i < 2 * pi; i += 2 * pi / 1000) {
            x.push_back(i);
            y.push_back(-cos(i) + 1);
        }
    }
    else {
        for (double i = 0; i < 10; i += 1.0 / 1000) {
            x.push_back(i);
            y.push_back(i * i * exp(-i));
        }
    }
    vector<double> orig = y;

    double height = 0.3;
    int location;
    int width = 5;
    int nums = 5;
    for (int i = 0; i < nums; i++) {
        //get the index of the spike
        int index = rand() % 1000;
        //get the height
        double lower = 0.0;
        double upper = height;
        std::uniform_real_distribution<double> unif(lower, upper);
        std::default_random_engine re;
        double h = unif(re);
        y[index] += h;
        //set the neightbors near the spike
        if (width > 1) {
            double left_range = y[index] - y[index - width / 2 - 1];
            double right_range = abs(y[index + width / 2 + 1] - y[index]);
            double l = 0.0;
            for (int i = index - width / 2; i < index; i++) {
                double u = left_range;
                //std::uniform_real_distribution<double> uni(lower, upper);
                //std::default_random_engine ran;
                //double change = uni(ran);
                double change = fRand(l, u);
                y[i] += change;
                l = change;
            }
            double u = right_range;
            for (int i = index + 1; i < index + width / 2 + 1; i++) {
                double l = 0.0;
                //std::uniform_real_distribution<double> uni(lower, upper);
                //std::default_random_engine ran;
                //double change = uni(ran);
                double change = fRand(l, u);
                y[i] += change;
                u = change;
            }
        }
    }
    vector<double> add_spike = y;

    spike_elim(x, y, 3, 1.2);
    double orig_sum = 0.0;
    double spike_sum = 0.0;
    double remove_sum = 0.0;
    for (double d : orig)
        orig_sum += d;
    for (double d : add_spike)
        spike_sum += d;
    for (double d : y)
        remove_sum += d;
    file1 << "x, y_orig, y_spike, y_remove" << endl;
    for (int i = 0; i < 1000; i++) {
        file1 << x[i] << ", ";
        file1 << orig[i] << ", ";
        file1 << add_spike[i] << ", ";
        file1 << y[i] << "," << endl;
    }
    return 0;
}

