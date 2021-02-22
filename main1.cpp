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
#include "File_Manager.hpp"
#include "FileHandler.hpp"
#include "smartPeakDetect.hpp"
#include "DataSmoothing.hpp"
#include "Levenberg_Marquardt.hpp"
#include "quick_half_max.hpp"
#include "data_input.hpp"
#include "remove_spike.hpp"
#include "Savitzky_Golay.hpp"
#include "background_remove.hpp"
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
//if unistd.h can't be included then comment out the else clause
//#else
//#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;

double func(const double input, const vector<double> params) {
    double T = 0.0;
    double I_t = 0.0;
    double k = .000086173303;
    double energy = params[0];
    double Tm = params[1] + 273.15;
    double dm = (2.0 * k * (Tm)) / energy;
    double Im = params[2];
    T = double(input + 273.15);
    I_t = Im * exp(1.0 + (energy / (k * T)) * ((T - Tm) / Tm) - ((T * T) / (Tm * Tm)) * exp((energy / (k * T)) * ((T - Tm) / Tm)) * (1.0 - ((2.0 * k * T) / energy)) - dm);
    if (isinf(I_t)) {
        cout << T << " " << Tm << " " << Im << " " << dm << endl;
    }
    return I_t;
}

void gd(const vector<double>& temp, const vector<double>& curve, vector<vector<double>>& peakParams, double& FOM, ofstream& f1, bool& one) {
    //cout << peakParams.size() << endl;
    //temperary vector to store peak data
    vector<vector<double>> temp_params = peakParams;
    //temperary vector to store accumulated fitted count
    vector<double> temp_output(curve.size(), 0.0);
    int curveSize = int(curve.size());
    int peakNum = int(peakParams.size());
    double k = .000086173303;
    double rate1 = 0.0000001;
    //double rate1 = 0.0000000047;
    double rate2 = 0.00002;
    //double rate2 = 0.000001;
    //double rate3 = 0.001;
    double rate3 = 0.00001;
    //double rate3 = 0.0000065;
    double current_FOM = FOM;
    int iteration = 0;
    bool check = true;
    while (iteration < 500 && check) {
        vector<vector<double>> update(peakNum, vector<double>(3, 0.0));
        vector<vector<double>> curves(temp_params.size(), vector<double>(curve.size()));
        vector<double> totalCurve(curve.size());
        for (int i = 0; i < int(curve.size()); ++i) {
            double out;
            double accum = 0.0;
            for (int x = 0; x < int(temp_params.size()); ++x) {
                out = quickFok(temp[i], temp_params[x]);
                curves[x][i] = out;
                accum += out;
            }
            totalCurve[i] = accum;
        }
        //use FWHM to find the left and right half max points
        for (int b = 0; b < peakNum; b++) {
            int TL = temp_params[b][3];
            int TR = temp_params[b][5];
            double energy = temp_params[b][0];
            double Tm = temp_params[b][1];
            double Im = temp_params[b][2];
            for (int index = TL; index < TR + 1; index++) {
            //for (int index = 0; index < int(curve.size()); index++) {
                double y = curve[index] - totalCurve[index] + curves[b][index];
                if (y < 0)
                    y = 0;
                double T = temp[index];
                double deriv_E = -2.0 * Im * (-(2.0 * k * T * T * T * exp((energy * (T - Tm)) / (Tm * k * T))) / (energy * energy * Tm * Tm) + (2.0 * Tm * k) /
                    (energy * energy) - (T * (T - Tm) * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm * Tm * k) + (T - Tm) /
                    (Tm * k * T)) * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) /
                    (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) /
                        (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0));
                double deriv_Tm = -2.0 * Im * ((2.0 * T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm * Tm) - (T * T * (1.0 - (2.0 * k * T) / energy)
                    * exp((energy * (T - Tm)) / (Tm * k * T)) * (-(energy * (T - Tm)) / (Tm * Tm * k * T) - energy / (Tm * k * T))) / (Tm * Tm) - (energy * (T - Tm)) / (Tm * Tm * k * T) -
                    energy / (Tm * k * T) - (2.0 * k) / energy) * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) -
                    (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) -
                        (2.0 * Tm * k) / energy + 1.0));
                double deriv_Im = -2.0 * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) /
                    (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) /
                    (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0));
                update[b][0] += rate1 * deriv_E;
                update[b][1] += rate2 * deriv_Tm;
                update[b][2] += rate3 * deriv_Im;
            }
        }
        //vector<double> change(peakParams.size(), -1000);
        //for (int i = 0; i < peakNum; i++) {
        //    cout << update[i][0] << "   " << update[i][1] << "   " << update[i][2] << endl;
        //    if (abs(update[i][0]) > change[0])
        //        change[0] = abs(update[i][0]);
        //    if (abs(update[i][1]) > change[1])
        //        change[1] = abs(update[i][1]);
        //    if (abs(update[i][2]) > change[2])
        //        change[2] = abs(update[i][2]);
        //    //if (abs(update[i][2]) > change[2])
        //    //    change[2] = update[i][2];
        //    if (abs(change[0]) < 0.00001 || abs(change[1]) < 0.0001 || abs(change[2]) < 0.0001) {
        //    //if (abs(change[0]) < 0.00001 || change[2] > 0) {
        //        check = false;
        //    }
        //}
        //cout << endl;
        //apply the change to the peak paramter
        for (int c = 0; c < peakNum; c++) {
            temp_params[c][0] -= update[c][0];
            temp_params[c][1] -= update[c][1];
            temp_params[c][2] -= update[c][2];
        }
        //for (int j = 0; j < peakNum; j++) {
        //    cout << update[j][0] << ", " << update[j][1] << ", " << update[j][2] << endl;
        //    cout << temp_params[j][0] << ", " << temp_params[j][1] << ", " << temp_params[j][2] << endl;
        //}
        double new_integral = 0.0;
        vector<double> new_output(curve.size(), 0.0);
        for (int d = 0; d < curveSize; d++) {
            double new_fit = 0.0;
            for (int e = 0; e < peakNum; e++) {
                new_fit += func(temp[d], temp_params[e]);
            }
            new_integral += new_fit;
            new_output[d] = new_fit;
        }
        double new_fom = 0.0;
        for (int f = 0; f < curveSize; ++f) {
            new_fom += abs(curve[f] - new_output[f]) / new_integral;
        }
        if (new_fom > current_FOM)
            check = false;
        else
            current_FOM = new_fom;
        cout << ".";
        cout.flush();
        iteration++;
    }
    if (iteration == 1)
        one = true;
    FOM = current_FOM;
    peakParams = temp_params;
    //cout << "GD ran " << iteration << " times" << endl;
    //cout << "the new FOM is: " << FOM;
    f1 << " new fom: " << FOM << " " << iteration << " times" << endl;
}

int gd_types(const vector<double>& temp, const vector<double>& curve, vector<vector<double>>& peakParams, double& FOM, int type) {
    vector<vector<double>> change_range;
    if (type == 100)
        change_range = { {0.1, 0.05, 0.05}, {0.15, 0.04, 0.05}, {0.16, 0.03, 0.05} };
    else if (type == 200)
        change_range = { {0.16, 0.03, 0.05}, {0.18, 0.02, 0.05}, {0.4, 0.02, 0.05}, {0.3, 0.02, 0.05}, {0.15, 0.03, 0.05}, {0.2, 0.013, 0.05}, {0.18, 0.02, 0.05}, {0.2, 0.02, 0.05}, {0.18, 0.02, 0.05} };
    else if (type == 300)
        change_range = { {0.3, 0.03, 0.05}, {0.23, 0.02, 0.05}, {0.4, 0.03, 0.05}, {0.37, 0.03, 0.05}, {0.24, 0.023, 0.05}, {0.2, 0.012, 0.05}, {0.3, 0.024, 0.05}, {0.24, 0.023, 0.05} };
    else if (type == 400)
        change_range = { {0.15, 0.014, 0.05}, {0.21, 0.015, 0.05}, {0.18, 0.013, 0.05} };
    else if (type == 700)
        change_range = { {0.15, 0.03, 0.05}, {0.15, 0.017, 0.05}, {0.14, 0.02, 0.05}, {0.06, 0.013, 0.05} };
    else
        change_range = { {0.12, 0.035, 0.05}, {0.17, 0.016, 0.05}, {0.24, 0.08, 0.05}, {0.14, 0.026, 0.05}, {0.155, 0.012, 0.05}, {0.15, 0.02, 0.05}, {0.25, 0.053, 0.05} };
    int curveSize = int(curve.size());
    int peakNum = int(peakParams.size());
    //find the half max index
    vector<double> temperature = temp;
    find_index(temperature, peakParams);
    //calculate original integral
    double orig_integral = 0.0;
    vector<double> orig_fit(curve.size());
    for (int d = 0; d < int(curve.size()); d++) {
        double fit = 0.0;
        for (int e = 0; e < curveSize; e++) {
            fit += func(temp[d], peakParams[e]);
        }
        orig_integral += fit;
        orig_fit[d] = fit;
    }
    double orig_fom = 0.0;
    for (int f = 0; f < int(curve.size()); ++f) {
        orig_fom += abs(curve[f] - orig_fit[f]) / orig_integral;
    }
    vector<vector<double>> temp_params = peakParams;
    double k = .000086173303;
    double rate1 = 0.0000001;
    double rate2 = 0.00002;
    double rate3 = 0.00001;
    double current_FOM = orig_fom;
    int iteration = 0;
    //bool matrix to check if should apply change
    vector<vector<bool>> process(curveSize, vector<bool>(3, true));
    //check fom
    bool check = true;
    //check if any still inside allowable range
    bool continue_process = true;
    while (iteration < 500 && check && continue_process) {
        vector<vector<double>> update(peakNum, vector<double>(3, 0.0));
        vector<vector<double>> curves(temp_params.size(), vector<double>(curve.size()));
        vector<double> totalCurve(curve.size());
        for (int i = 0; i < int(curve.size()); ++i) {
            double out;
            double accum = 0.0;
            for (int x = 0; x < int(temp_params.size()); ++x) {
                out = quickFok(temp[i], temp_params[x]);
                curves[x][i] = out;
                accum += out;
            }
            totalCurve[i] = accum;
        }
        //use FWHM to find the left and right half max points
        for (int b = 0; b < peakNum; b++) {
            int TL = temp_params[b][3];
            int TR = temp_params[b][5];
            double energy = temp_params[b][0];
            double Tm = temp_params[b][1];
            double Im = temp_params[b][2];
            for (int index = TL; index < TR + 1; index++) {
                double y = curve[index] - totalCurve[index] + curves[b][index];
                if (y < 0)
                    y = 0;
                double T = temp[index];
                double deriv_E = -2.0 * Im * (-(2.0 * k * T * T * T * exp((energy * (T - Tm)) / (Tm * k * T))) / (energy * energy * Tm * Tm) + (2.0 * Tm * k) /
                    (energy * energy) - (T * (T - Tm) * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm * Tm * k) + (T - Tm) /
                    (Tm * k * T)) * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) /
                    (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) /
                        (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0));
                double deriv_Tm = -2.0 * Im * ((2.0 * T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm * Tm) - (T * T * (1.0 - (2.0 * k * T) / energy)
                    * exp((energy * (T - Tm)) / (Tm * k * T)) * (-(energy * (T - Tm)) / (Tm * Tm * k * T) - energy / (Tm * k * T))) / (Tm * Tm) - (energy * (T - Tm)) / (Tm * Tm * k * T) -
                    energy / (Tm * k * T) - (2.0 * k) / energy) * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) -
                    (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) -
                        (2.0 * Tm * k) / energy + 1.0));
                double deriv_Im = -2.0 * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) /
                    (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) /
                    (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0));
                update[b][0] += rate1 * deriv_E;
                update[b][1] += rate2 * deriv_Tm;
                update[b][2] += rate3 * deriv_Im;
            }
        }
        //apply changes and check if out range
        for (int c = 0; c < peakNum; c++) {
            if (process[c][0]) {
                temp_params[c][0] -= update[c][0];
                if ((abs(temp_params[c][0] - peakParams[c][0]) / peakParams[c][0]) > change_range[c][0]) {
                    //revert
                    temp_params[c][0] += update[c][0];
                    //stop process it
                    process[c][0] = false;
                }
            }
            if (process[c][1]) {
                temp_params[c][1] -= update[c][1];
                if ((abs(temp_params[c][1] - peakParams[c][1]) / peakParams[c][1]) > change_range[c][1]) {
                    //revert
                    temp_params[c][1] += update[c][1];
                    //stop process it
                    process[c][1] = false;
                }
            }
            if (process[c][2]) {
                temp_params[c][2] -= update[c][2];
                if ((abs(temp_params[c][2] - peakParams[c][2]) / peakParams[c][2]) > change_range[c][2]) {
                    //revert
                    temp_params[c][2] += update[c][2];
                    //stop process it
                    process[c][2] = false;
                }
            }
        }
        //check if all process is false
        continue_process = false;
        for (vector<bool>& b : process) {
            for (int i = 0; i < 3; i++) {
                if (b[i])
                    continue_process = true;
            }
        }
        //calculate new fom
        double new_integral = 0.0;
        vector<double> new_output(curve.size(), 0.0);
        for (int d = 0; d < curveSize; d++) {
            double new_fit = 0.0;
            for (int e = 0; e < peakNum; e++) {
                new_fit += func(temp[d], temp_params[e]);
            }
            new_integral += new_fit;
            new_output[d] = new_fit;
        }
        double new_fom = 0.0;
        for (int f = 0; f < curveSize; ++f) {
            new_fom += abs(curve[f] - new_output[f]) / new_integral;
        }
        if (new_fom > current_FOM)
            check = false;
        else
            current_FOM = new_fom;
        cout << ".";
        cout.flush();
        iteration++;
    }
    FOM = current_FOM;
    peakParams = temp_params;
    return iteration;
}

int main(int argc, char* argv[]) {
    //enable quick mode to run output with machine generated peak detections
    bool output_mode = false;
    if (argc > 1)
        output_mode = (*(argv[1]) == 'q') ? true : false;
    //string to store the input directory path
    string dir;
    //functional string that will be used in the while loop to read in path
    string start = "n";
    //vector to store the processed files, will be used in the statistic function call
    vector<string> filenames;
    //vector of strings that will store all the csv file names in the "dir" directory
    vector<string> files;
    string output_dir = "";
    //sampling rate default to be 0.1
    int time = 0.1;

    //read in the path that contains all the input data
    while (start == "n" || start == "N") {
        cout << "Please enter the full path to directory containing csv formatted emission spectra:" << endl;
        cin >> dir;
        //cast dir in the right format to be opened
        if (dir.back() == '/')
            dir.pop_back();

        //HANDLE FILES
        //call the handle_dir function from FileHandler.cpp, store all csv files in dir to the files vector
        files = handle_dir(dir, output_dir);
        if (!output_mode) {
            cout << "Is this correct and would you like to start processing (y/n)?" << endl;
            cin >> start;
        }
        else
            start = "y";
    }
    //2d vector to store processed data of each input file
    vector<vector<double>> stats(files.size(), vector<double>(0, 0.0));
    int count = 0;

    vector<vector<double>> peakParams;
    vector<vector<vector<double>>> all_peakParam;
    vector<int> samp_rate;
    enum Mode { ALL, EACH, NONE };
    Mode m = Mode::NONE;
    if (!output_mode) {
        //ASK USER PEAK INPUT
        string choice;
        cout << "Do you want to manually input peak locations and heights (y/n)?" << endl;
        cin >> choice;
        if (choice == "y") {
            string repeat;
            cout << "Do you want the same input for all files or different input for each file (all/each)?" << endl;
            cin >> repeat;
            //ask user to input single set of peaks data for all files
            if (repeat == "all") {
                m = Mode::ALL;
                cout << "Please type in data in the format: tmeperature,count,activation energy, press enter for each peak." << endl;
                cout << "Type done when you are finished." << endl;
                peakParams = input_data();
                if (peakParams.empty()) {
                    m = Mode::NONE;
                    cout << "Empty input, switching to automatic peak identification." << endl;
                }
                cout << "Please give the sampling rate for all data.";
                cin >> time;
            }
            else if (repeat == "each") {
                m = Mode::EACH;
                cout << "For each file, please type in peak data in the format: tmeperature,count,activation energy, press enter for each peak." << endl;
                cout << "Type done when you are finished." << endl;
                cout << "Then, please type in the sampling rate." << endl;
                auto i = files.begin();
                for (; i != files.end(); ++i) {
                    string filename = i->substr((i->find_last_of("/\\")) + 1);
                    cout << "Enter peak data for: " << filename << endl;
                    vector<vector<double>> param = input_data();
                    all_peakParam.push_back(param);
                    cout << "Enter sampling rate." << endl;
                    int t;
                    cin >> t;
                    samp_rate.push_back(t);
                }
                if (all_peakParam.empty()) {
                    m = Mode::NONE;
                    cout << "Empty input, switching to automatic peak identification." << endl;
                }
            }
            else {
                cout << "Invalid command!" << endl;
            }
        }
    }
    ofstream file1;
    string path = "C:/Users/jack0/Desktop/experiment.txt";
    file1.open(path);
    for (int i = 0; i < static_cast<int>(files.size()); i++) {
        //count to see which file is being processed
        int num = 0;
        //make a coopy of peakParam and then customize the copy accoding to the data read in
        vector<vector<double>> peak;
        if (m == Mode::ALL)
            peak = peakParams;
        if (m == Mode::EACH)
            peak = all_peakParam[num];

        //erase the previous temp.csv file to read in new data
        if (files[i].find("temp.csv") != string::npos) {
            files.erase(files.begin() + i);
            continue;
        }
        cout << "----------------------------" << endl << "Processing: ";
        string filename = files[i].substr((files[i].find_last_of("/\\")) + 1);
        cout << filename << " (" << count + 1 << " of " << files.size() << ")" << endl << "Reading in File  .";
        cout.flush();

        //FILE_MANAGER created
        //create a fileManager object that takes in the i/csv path
        File_Manager fileManager = *new File_Manager(files[i], output_dir);
        cout << ".";
        cout.flush();
        //create a pair of two vector data which has first to be temperature data and second to be count data
        pair<vector<double>, vector<double>> data = fileManager.read();
        int length = static_cast<int>(data.second.size());
        int window_size = length * 0.05;
        if (window_size % 2 == 0) {
            window_size += 1;
        }

        //Smooth temperature raw data with Savitzky-Golay
        int window = length * 0.05;
        if (window % 2 == 0) {
            window += 1;
        }
        vector<double> temp_orig1 = data.first;
        vector<double> temp_orig2 = temp_orig1;
        //SG_smooth(temp_orig1, window, 4);
        //SG_smooth(temp_orig2, window, 5);
        SG_smooth(temp_orig2, window, 3);
        //for (int i = 0; i < static_cast<int>(temp_orig1.size()); i++) {
        //    data.first[i] = (temp_orig1[i] + temp_orig2[i]) / 2;
        //}
        for (int i = 0; i < static_cast<int>(temp_orig1.size()); i++) {
            data.first[i] = temp_orig2[i];
        }

        vector<double> orig_count = data.second;
        //REMOVE_SPIKE call
        spike_elim(data.first, data.second, 3, 1.2);

        //copy two times the count data and run Savitzky-Golay with order 4 and 5, then take the average
        vector<double> orig_count1 = data.second;
        vector<double> orig_count2 = orig_count1;
        SG_smooth(orig_count1, window_size, 4);
        SG_smooth(orig_count2, window_size, 5);
        for (int i = 0; i < static_cast<int>(orig_count1.size()); i++) {
            data.second[i] = (orig_count1[i] + orig_count2[i]) / 2;
        }

        double max_intensity = *max_element(data.second.begin(), data.second.end());
        //background_substraction
        //vector<double> t = remove_back(data.first, data.second);
        vector<double> temp = data.second;

        //calculate the curve area by adding the count data
        const double curveArea = accumulate(data.second.begin(), data.second.end(), 0.0);
        //if the curve area is less 2000 then it's not enough for further analysis
        if (curveArea < 2000) {
            files.erase(files.begin() + i);
            i--;
            cout << "file not considered" << endl;
            //remove((dir + "/temp.csv").c_str());
            remove((output_dir + "/temp.csv").c_str());
            continue;
        }

        //testing gradient descent on TLD 100
        //peakparams: activation energy, maxTemp, maxIntensity, TL, TM, TR
        //vector<vector<double>> peakParams;
        //findPeaks(data.first, data.second, peakParams, output_dir);
        //if (int(peakParams.size()) != 4) {
        //    files.erase(files.begin() + i);
        //    i--;
        //    cout << "file not considered" << endl;
        //    //remove((dir + "/temp.csv").c_str());
        //    remove((output_dir + "/temp.csv").c_str());
        //    continue;
        //}
        //cout << "energy, temp, count, TL, TM, TR" << endl;
        //for (int i = 0; i < int(peakParams.size()); i++)
        //    cout << peakParams[i][0] << " " << peakParams[i][1] << " " << peakParams[i][2] << " " << peakParams[i][3] << " " << peakParams[i][4] << " " << peakParams[i][5] << endl;
        //vector<vector<double>> curve;
        //for (int i = 0; i < int(peakParams.size()); ++i) {
        //    curve.push_back(vector<double>(data.first.size(), 0.0));
        //}
        //calculate every temperature's FOK data in each peak fit, accumulate peak areas for each peak in
        //peak_areas and accumulate same temperature's FOK values in all fits to sum

        //TLD 300
        vector<vector<double>> peak_param;
        //peak_param.push_back({1.45, 93.85, 0.51});
        //peak_param.push_back({ 1.17, 111.85, 1.2});
        //peak_param.push_back({ 1.3, 132.85, 1.8});
        //peak_param.push_back({ 1.09, 160.85, 1.25});
        //peak_param.push_back({ 1.31, 187.85, 7.1});
        //peak_param.push_back({ 1.26, 215.85, 0.5});
        //peak_param.push_back({ 1.43, 264.85, 1.85});
        //peak_param.push_back({ 1.31, 288.85, 1.75});
        //TLD 400
        //peak_param.push_back({ 1.42, 296.85, 4});
        //peak_param.push_back({ 1.32, 323.85, 6.2});
        //peak_param.push_back({ 1.36, 349.85, 4.8});
        //TLD 900
        //peak_param.push_back({ 1.02, 122.85, 0.38});
        //peak_param.push_back({ 0.96, 149.85, 0.40});
        //peak_param.push_back({ 1.05, 161.85, 0.65});
        //peak_param.push_back({ 1.40, 188.85, 0.55});
        //peak_param.push_back({ 1.42, 205.85, 0.25});
        //peak_param.push_back({ 0.96, 249.85, 0.65});
        //peak_param.push_back({ 0.88, 272.85, 3.25});
        //vector<double> total_curve(data.first.size());

        //calculate original integral
        //double orig_integral = 0.0;
        //vector<double> orig_fit(data.second.size());
        //for (int d = 0; d < int(data.second.size()); d++) {
        //    double fit = 0.0;
        //    for (int e = 0; e < int(peak_param.size()); e++) {
        //        fit += func(data.first[d], peak_param[e]);
        //    }
        //    orig_integral += fit;
        //    orig_fit[d] = fit;
        //}
        //double orig_fom = 0.0;
        //for (int f = 0; f < int(data.second.size()); ++f) {
        //    orig_fom += abs(data.second[f] - orig_fit[f]) / orig_integral;
        //}

        vector<vector<double>> peak_100;
        double fom_100 = -1;
        int iteration_100 = gd_types(data.first, data.second, peak_100, fom_100, 100);
        vector<vector<double>> peak_200;
        double fom_200 = -1;
        int iteration_200 = gd_types(data.first, data.second, peak_200, fom_200, 200);
        vector<vector<double>> peak_300;
        double fom_300 = -1;
        int iteration_300 = gd_types(data.first, data.second, peak_300, fom_300, 300);
        vector<vector<double>> peak_400;
        double fom_400 = -1;
        int iteration_400 = gd_types(data.first, data.second, peak_400, fom_400, 400);
        vector<vector<double>> peak_700;
        double fom_700 = -1;
        int iteration_700 = gd_types(data.first, data.second, peak_700, fom_700, 700);
        vector<vector<double>> peak_900;
        double fom_900 = -1;
        int iteration_900 = gd_types(data.first, data.second, peak_900, fom_900, 900);
        vector<double> fom_set{ fom_100 , fom_200, fom_300, fom_400, fom_700, fom_900 };
        int min_fom = min_element(fom_set.begin(), fom_set.end()) - fom_set.begin();
        if (min_fom == 0)
            peak_param = peak_100;
        else if (min_fom == 1)
            peak_param = peak_200;
        else if(min_fom == 2)
            peak_param = peak_300;
        else if(min_fom == 3)
            peak_param = peak_400;
        else if (min_fom == 4)
            peak_param = peak_700;
        else
            peak_param = peak_900;


        //double area = 0.0;
        //for (int i = 0; i < int(data.first.size()); ++i) {
        //    for (int x = 0; x < int(peak_param.size()); ++x) {
        //        area += quickFok(data.first[i], peak_param[x]);
        //    }
        //}
        //int cons = 1;
        //int iteration = 0;
        //vector<vector<double>> temp_param = peak_param;
        //while (area < 0.97 * curveArea || area > 1.03 * curveArea) {
        //    if (iteration > 1000)
        //        break;
        //    cons = 20;
        //    cons *= (1 - (area / curveArea));
        //    for (auto& v : temp_param)
        //        v[2] += cons;
        //    area = 0.0;
        //    for (int i = 0; i < int(data.first.size()); ++i) {
        //        for (int x = 0; x < int(peak_param.size()); ++x) {
        //            area += quickFok(data.first[i], temp_param[x]);
        //        }
        //    }
        //    iteration++;
        //}
        //double new_integral = 0.0;
        //vector<double> new_fit(data.second.size());
        //for (int d = 0; d < int(data.second.size()); d++) {
        //    double fit = 0.0;
        //    for (int e = 0; e < int(temp_param.size()); e++) {
        //        fit += func(data.first[d], temp_param[e]);
        //    }
        //    new_integral += fit;
        //    new_fit[d] = fit;
        //}
        //double new_fom = 0.0;
        //for (int f = 0; f < int(data.second.size()); ++f) {
        //    new_fom += abs(data.second[f] - new_fit[f]) / new_integral;
        //}
        //double cur_fom = 0.0;
        //for (int f = 0; f < int(total_curve.size()); ++f) {
        //    cur_fom += abs(data.second[f] - total_curve[f]) / total;
        //}
        //double final_cons = temp_param[0][2] / peak_param[0][2];
        //file1 << filename << " ";
        //file1 << "orig_fom: " << orig_fom << " new_fom: " << new_fom << " ratio: " << final_cons << " iterations: " << iteration << " max_intensity: " << max_intensity << endl;
        //vector<vector<double>> GDParams = peakParams;
        //vector<vector<double>> GDcurve;
        //double fom = 1;
        //vector<vector<double>> GDParams = peakParams;
        //bool check = false;
        //gd(data.first, data.second, GDParams, cur_fom, file1, check);
        
        //gd(data.first, data.second, GDParams, fom);
        //for (int i = 0; i < int(GDParams.size()); ++i) {
        //    GDcurve.push_back(vector<double>(data.first.size(), 0.0));
        //}
        //for (int i = 0; i < int(data.first.size()); ++i) {
        //    double output = 0.0;
        //    for (int x = 0; x < int(GDParams.size()); ++x) {
        //        double out = quickFok(data.first[i], GDParams[x]);
        //        GDcurve[x][i] = out;
        //    }
        //}
        ////cout << oldParams[0][0] << " " << oldParams[0][1] << " " << oldParams[0][2] << endl;
        //cout << peakParams[0][0] << " " << peakParams[0][1] << " " << peakParams[0][2] << endl;

        
        //calculate the noise ratio
        //file1 << filename << " ";
        //double sum_smooth = 0.0;
        //double sum_diff = 0.0;
        //for (int i = 0; i < static_cast<int>(orig_count.size()); i++) {
        //    sum_smooth += data.second[i];
        //    sum_diff += abs(data.second[i] - orig_count[i]);
        //}
        //file1 << sum_smooth / sum_diff << endl;
        //string output = files[i] + "_output.csv";
        //file.open(output);
        //string output = files[i].substr(files[i].find("R"));
        //string path = output_dir + "/" + output;
        //string path = output_dir + "/" + filename;
        //if (check) {
        //    ofstream file2;
        //    string path = output_dir + "/" + filename;
        //    file2.open(path);
        //    file2 << "temp, after removal, first, sec, third, forth, GDfirst, GDsec, GDthird, GDforth";
        //    file2 << ",\n";
        //    for (int i = 0; i < int(data.first.size()); i++) {
        //        file2 << data.first[i] << ",";
        //        //file2 << orig_count[i] << ",";
        //        file2 << temp[i] << ",";
        //        for (int j = 0; j < int(curve.size()); j++) {
        //            file2 << curve[j][i] << ",";
        //        }
        //        for (int j = 0; j < int(GDcurve.size()) - 1; j++) {
        //            file2 << GDcurve[j][i] << ",";
        //        }
        //        file2 << GDcurve[curve.size() - 1][i];
        //        file2 << ",\n";
        //        
        //    }
        //    file2.close();
        //}
        
        //file2 << "temp, orig_count, new_count, subtracted_count, deriv";
        //file2 << ",\n";

        

        //file2.setf(ios_base::fixed);
        //file2 << setprecision(5);
        //for (int i = 0; i < int(orig_count.size()); ++i) {
        //for (int i = 0; i < int(data.second.size()); ++i) {
        //    file2 << data.first[i] << ",";
        //    file2 << orig_count[i] << ", ";
        //    file2 << temp[i] << ", ";
        //    file2 << data.second[i] << ", ";
        //    //file2 << t[i];
        //    file2 << ",\n";
        //}

        
        
        count++;
    }
    //file1.close();
    return 0;
}

