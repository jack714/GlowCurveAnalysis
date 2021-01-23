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
#include "Savitzky_Golay.hpp"
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

void generate_spike(vector<double>& y, vector<double>& x, vector<int> locations, double height, int width) {
    for (int i : locations) {
        //get the index of the spike
        int index = i;
        //get the height
        //double lower = 0.0;
        //double upper = height;
        y[index] += height;
        //set the neightbors near the spike
        if (width > 1) {
            //double left_range = y[index] - y[index - width / 2 - 1];
            //double right_range = abs(y[index + width / 2 + 1] - y[index]);
            //double l = 0.0;
            //for (int j = index - width / 2; j < index; j++) {
            //    double u = left_range;
            //    double change = fRand(l, u);
            //    y[j] += change;
            //    l = change;
            //}
            //double u = right_range;
            //for (int k = index + 1; k < index + (width / 2) + 1; k++) {
            //    double l = 0.0;
            //    double change = fRand(l, u);
            //    y[k] += change;
            //    u = change;
            //}
            int left = index - width / 2 - 1;
            int right = index + width / 2 + 1;
            double left_slope = (y[index] - y[left]) / (x[index] - x[left]);
            double left_k = y[index] - x[left] * left_slope;
            double right_slope = (y[right] - y[index]) / (x[right] - x[index]);
            double right_k = y[index] - right_slope * x[index];
            for (int i = left + 1; i < index; i++) {
                y[i] = x[i] * left_slope + left_k;
            }
            for (int j = index + 1; j < right; j++) {
                y[j] = x[j] * right_slope + right_k;
            }
        }
    }
}

int main() {
    ofstream file1;
    string path = "C:/Users/jack0/Desktop/output1.csv";
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
    //int width = 5;
    file1 << "Orignal, Add spike, SG, spike_remove, " << endl;
    for (int i = 0; i < 10000; i += 1) {
        int width = rand() % 13;
        vector<int> locations;
        //locations.push_back(500);
        double d = 0.0;
        while (d < 5.0) {
            locations.push_back(floor(fRand(100, 900)));
            d = fRand(0, 7.0);
        }

        //double height = 4.0;
        double height = fRand(0, 8.0);
        
    
        //ofstream file2;
        //string p = "C:/Users/jack0/Desktop/" + to_string(h) + ".csv";
        //file2.open(p);
        //height += d;
        //width = w;
        //height = h;
        //vector<int> loc;
        //loc.push_back(i);
        vector<double> temp = orig;
        generate_spike(temp, x, locations, height, width);
        vector<double> spike = temp;

        vector<double> orig_copy_1 = temp;
        vector<double> orig_copy_2 = temp;
        vector<double> sg_copy = temp;
        int window_size = 201;
        SG_smooth(orig_copy_1, window_size, 4);
        SG_smooth(orig_copy_2, window_size, 5);
        for (int i = 0; i < 1000; i++) {
            sg_copy[i] = (orig_copy_1[i] + orig_copy_2[i]) / 2;
        }

        vector<double> add_spike = temp;
        spike_elim(x, temp, 3, 1.2);
        double orig_sum = 0.0;
        double spike_sum = 0.0;
        double remove_sum = 0.0;
        double sg_sum = 0.0;
        for (double d : orig)
            orig_sum += d;
        for (double d : add_spike)
            spike_sum += d;
        for (double d : temp)
            remove_sum += d;
        for (double d : sg_copy)
            sg_sum += d;
        
        file1 << orig_sum << ", ";
        file1 << spike_sum << ", ";
        file1 << sg_sum << ", ";
        file1 << remove_sum << ",\n";

        //file2 << "temp, orig, spike, sg, remove" << endl;
        //for (int i = 0; i < 1000; i++) {
        //    file2 << x[i] << ", ";
        //    file2 << orig[i] << ", ";
        //    file2 << spike[i] << ", ";
        //    file2 << sg_copy[i] << ", ";
        //    file2 << temp[i] << ",\n";
        //}

    }
    //vector<double> orig_copy_1 = y;
    //vector<double> orig_copy_2 = y;
    //vector<double> sg_copy = y;
    //int window_size = 201;
    //SG_smooth(orig_copy_1, window_size, 4);
    //SG_smooth(orig_copy_2, window_size, 5);
    //for (int i = 0; i < 1000; i++) {
    //    sg_copy[i] = (orig_copy_1[i] + orig_copy_2[i]) / 2;
    //}
    //vector<double> add_spike = y;
    //
    //spike_elim(x, y, 3, 1.2);
    //
    //double orig_sum = 0.0;
    //double spike_sum = 0.0;
    //double remove_sum = 0.0;
    //double sg_sum = 0.0;
    //for (double d : orig)
    //    orig_sum += d;
    //for (double d : add_spike)
    //    spike_sum += d;
    //for (double d : y)
    //    remove_sum += d;
    //for (double d : sg_copy)
    //    sg_sum += d;
    //cout << "Original: " << orig_sum << endl;
    //cout << "Add spike: " << spike_sum << endl;
    //cout << "SG: " << sg_sum << endl;
    //cout << "spike remove: " << remove_sum << endl;

    //file1 << "x, y_orig, y_spike, y_SG, y_remove" << endl;
    //for (int i = 0; i < 1000; i++) {
    //    file1 << x[i] << ", ";
    //    file1 << orig[i] << ", ";
    //    file1 << add_spike[i] << ", ";
    //    file1 << sg_copy[i] << ", ";
    //    file1 << y[i] << "," << endl;
    //}
    return 0;
}

