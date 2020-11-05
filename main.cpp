//
//  main.cpp
//  GlowCurveAnalsys
//
//  Initially created by jeremy hepker on 1/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//

#include <iostream>
#include <getopt.h>
#include <string>
#include <locale>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include "File_Manager.hpp"
#include "FileHandler.hpp"
#include "smartPeakDetect.hpp"
#include "DataSmoothing.hpp"
#include "Levenberg-Marquardt.hpp"
#include "quick_half_max.hpp"
#include "data_input.hpp"
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
    bool output_mode = (*(argv[1]) == 'q') ? true : false;
    cout << *(argv[1]);
    //string to store the input directory path
    string dir;
    //functional string that will be used in the while loop to read in path
    string start = "n";
    //vector to store the processed files, will be used in the statistic function call
    vector<string> filenames;
    //vector of strings that will store all the csv file names in the "dir" directory
    vector<string> files;
    string output_dir = "";
    
    //read in the path that contains all the input data
    while(start == "n" || start =="N"){
        cout<< "Please enter the full path to directory containing csv formatted emission spectra:" << endl;
        cin >> dir;
        //cast dir in the right format to be opened
        if(dir.back() == '/') 
            dir.pop_back();
        
        //HANDLE FILES
        //call the handle_dir function from FileHandler.cpp, store all csv files in dir to the files vector
        files = handle_dir(dir, output_dir);
        cout << "Is this correct and would you like to start processing (y/n)?" << endl;
        cin >> start;
    }
    //2d vector to store processed data of each input file
    vector<vector<double>> stats(files.size(), vector<double>(0,0.0));
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
            }
            else {
                cout << "Invalid command!" << endl;
            }
        }
    }
    auto i = files.begin();
    for (; i != files.end(); ++i) {
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
        if (i->find("temp.csv") != string::npos) {
            files.erase(i);
            continue;
        }
        cout << "----------------------------" << endl << "Processing: ";
        string filename = i->substr((i->find_last_of("/\\")) + 1);
        cout << filename << " (" << count + 1 << " of " << files.size() << ")" << endl << "Reading in File  .";
        cout.flush();

        //FILE_MANAGER created
        //create a fileManager object that takes in the i/csv path
        File_Manager fileManager = *new File_Manager(*i);
        cout << ".";
        cout.flush();
        //create a pair of two vector data which has first to be temperature data and second to be count data
        pair<vector<double>, vector<double>> data = fileManager.read();

        //DATA_SMOOTHING call
        //use dataSmooth from dataSmoothing.cpp to process raw data
        for (int i = 0; i < 5; ++i)
            dataSmooth(data.first, data.second);
        //calculate the curve area by adding the count data
        const double curveArea = accumulate(data.second.begin(), data.second.end(), 0.0);
        //if the curve area is less 2000 then it's not enough for further analysis
        if (curveArea < 2000) {
            files.erase(i);
            continue;
        }
        //remove the temp.csv created in fileManager since already read them in data
        remove((dir + "/temp.csv").c_str());

        //if it's the quick output mode or user didn't input data, let the program detect peaks
        if (output_mode || m == Mode::NONE) {
            //SMART_PEAK_DETECTION
            cout << "." << endl << "Finding Peaks  ..";
            cout.flush();
            //call findPeaks from smartPeakDetect.cpp
            firstDir = findPeaks(data.first, data.second, peak, output_dir);
            cout << ".";
            cout.flush();
            stats[count].push_back(fileManager.barcode());
            //need cast to cstr for remove in stdio.h
            remove((dir + "/temp.csv").c_str());
            cout.flush();
        }

        //populate peakParams with full width half max indexes
        //find_index(data.first, peakParams);

        //LEVENBERG-MARQUART call
        cout << endl << "Deconvoluting Glow Peak  .";
        cout.flush();
        //calling Levenberg-Marquardt.cpp and create a First_Order_Kinetics object
        First_Order_Kinetics FOK_Model = *new First_Order_Kinetics(data, peak);
        //calculate FOM and the area under each curve fit
        stats[count].push_back(FOK_Model.glow_curve());
        stats[count].push_back(FOK_Model.area());
        //push fileManager's heating rate to back of stats
        stats[count].push_back(fileManager.temp_rate(filename));
        //copy FOK_Model's curve_area vector which contains area under curve for each peak to peak_integral
        vector<double> peak_integral = FOK_Model.return_curve_areas();
        //store the area under curve for each peak fit to stats
        for (auto i = peak_integral.begin(); i != peak_integral.end(); ++i) {
            stats[count].push_back(*i);
        }

        //OUTPUT to files
        //assign returnedPeak to be 2d vector that contains FOK data for each peak fit
        vector<vector<double>> returnedPeaks = FOK_Model.return_glow_curve();
        filenames.push_back(filename);
        filename = output_dir + "/" + filename;
        //create output csv file and write temperature, count, fitted count data
        fileManager.write(returnedPeaks, filename);
        cout << "----------------------------" << endl;
        count++;
        //when all files are read, output statistic file for an overview of all fittings
        if (count == int(files.size()))
            fileManager.statistics(stats, filenames, output_dir);
        num++;
    }
    return 0;
}

