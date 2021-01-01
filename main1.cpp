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
#include "background_subtraction.hpp"
#include "remove_spike.hpp"
#include "Savitzsky_Golay.hpp"
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
//if unistd.h can't be included then comment out the else clause
//#else
//#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;

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
    //string path = "C:/Users/wenji/Desktop/LARGE_GCA_TEST_SUITE/large_output.txt";
    //file1.open(path);
    for (int i = 0; i < static_cast<int>(files.size()); i++) {
        //count to see which file is being processed
        int num = 0;
        vector<double> firstDir;
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
        vector<double> orig_count = data.second;
        //REMOVE_SPIKE call
        spike_elim(data.first, data.second, 3, 1.2);
        //DATA_SMOOTHING call
        //use dataSmooth from dataSmoothing.cpp to process raw data
        //for (int j = 0; j < 5; ++j)
        //    dataSmooth(data.first, data.second);
        //dataSmooth(data.first, data.second);

        //copy two times the count data and run Savitzsky-Golay with order 4 and 5, then take the average
        vector<double> orig_count1 = data.second;
        vector<double> orig_count2 = orig_count1;
        SG_smooth(orig_count1, window_size, 4);
        SG_smooth(orig_count2, window_size, 5);
        for (int i = 0; i < static_cast<int>(orig_count1.size()); i++) {
            data.second[i] = (orig_count1[i] + orig_count2[i]) / 2;
        }
        //calculate the curve area by adding the count data
        const double curveArea = accumulate(data.second.begin(), data.second.end(), 0.0);
        //if the curve area is less 2000 then it's not enough for further analysis
        if (curveArea < 2000) {
            files.erase(files.begin() + i);
            i--;
            cout << endl;
            //remove((dir + "/temp.csv").c_str());
            remove((output_dir + "/temp.csv").c_str());
            continue;
        }
        ofstream file2;
        //calculate the noise ratio
        file1 << filename << " ";
        double sum_smooth = 0.0;
        double sum_diff = 0.0;
        for (int i = 0; i < static_cast<int>(orig_count.size()); i++) {
            sum_smooth += data.second[i];
            sum_diff += abs(data.second[i] - orig_count[i]);
        }
        file1 << sum_smooth / sum_diff << endl;
        //string output = files[i] + "_output.csv";
        //file.open(output);
        //string output = files[i].substr(files[i].find("R"));
        //string path = output_dir + "/" + output;
        string path = output_dir + "/" + filename;
        file2.open(path);
        file2 << "temp, orig_count, new_count";
        file2 << ",\n";
        file2.setf(ios_base::fixed);
        file2 << setprecision(5);
        //for (int i = 0; i < int(orig_count.size()); ++i) {
        for (int i = 0; i < int(data.second.size()); ++i) {
            file2 << data.first[i] << ",";
            file2 << orig_count[i] << ", ";
            file2 << data.second[i];
            file2 << ",\n";
        }
        file2.close();
        count++;
    }
    file1.close();
    return 0;
}

