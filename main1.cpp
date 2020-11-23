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
    //record the time elapse between data points
    double time = 0.0;

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
            }
            else if (repeat == "each") {
                m = Mode::EACH;
                cout << "For each file, please type in data in the format: tmeperature,count,activation energy, press enter for each peak." << endl;
                cout << "Type done when you are finished." << endl;
                auto i = files.begin();
                for (; i != files.end(); ++i) {
                    string filename = i->substr((i->find_last_of("/\\")) + 1);
                    cout << "Enter data for: " << filename << endl;
                    vector<vector<double>> param = input_data();
                    all_peakParam.push_back(param);
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
        pair<vector<double>, vector<double>> data = fileManager.read(time);
        vector<double> orig_count = data.second;
        //REMOVE_SPIKE call
        spike_elim(data.first, data.second, 3, 1.2);
        //DATA_SMOOTHING call
        //use dataSmooth from dataSmoothing.cpp to process raw data
        for (int j = 0; j < 5; ++j)
            dataSmooth(data.first, data.second);
        ofstream file;
        string output = files[i] + "_output2.csv";
        file.open(output);
        file << "temp, orig_count, new_count";
        file << ",\n";
        file.setf(ios_base::fixed);
        file << setprecision(5);
        //for (int i = 0; i < int(orig_count.size()); ++i) {
        for (int i = 0; i < int(data.second.size()); ++i) {
            file << data.first[i] << ",";
            file << orig_count[i] << ", ";
            file << data.second[i];
            file << ",\n";
        }
        file.close();
    }
    return 0;
}

