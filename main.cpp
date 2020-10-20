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
#include "Levenberg–Marquardt.hpp"
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
//if unistd.h can't be included then comment out the else clause
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;

int main() {
    //string to store the input directory path
    string dir;
    //functional string that will be used in the while loop to read in path
    string start ="n";
    //vector to store the processed files, will be used in the statistic function call
    vector<string> filenames;
    //vector of strings that will store all the csv files in the "dir" directory
    vector<string> files;
    string output_dir = "";
    
    //read in the path that contains all the input data
    while(start =="n" || start =="N"){
        cout<< "Please enter the full path to directory containing csv formatted emission spectra:" << endl;
        cin >> dir;
        //cast dir in the right format to be opened
        if(dir.back() == '/') 
            dir.pop_back();
        //call the handle_dir function from FileHandler cpp
        //store all csv files in dir to the files vector
        files = handle_dir(dir, output_dir);
        cout << "Is this correct and would you like to start processing (y/n)?"<<endl;
        cin >> start;
    }
    //2d vector to store processed data of each input file
    vector<vector<double>> stats(files.size(), vector<double>(0,0.0));
    int count = 0;
    auto i = files.begin();
    //iterate all the csv files in files and process them
    for(; i != files.end(); ++i){
        vector<vector<double>> peakParams;
        double integral= 0;
        //erase the previous temp.csv file to read in new data
        if(i->find("temp.csv") != string::npos){
            files.erase(i);
            continue;
        }
        cout << "----------------------------" << endl << "Processing: ";
        string filename = i->substr((i->find_last_of("/\\"))+1);
        cout << filename << " (" << count+1 << " of " << files.size() << ")" << endl << "Reading in File  .";
        cout.flush();
        //create a fileManager object that takes in the i/csv path
        File_Manager fileManager = *new File_Manager(*i);
        cout<<".";
        cout.flush();
        pair<vector<double>, vector<double>> data = fileManager.read();
        const double curveArea = accumulate(data.second.begin(), data.second.end(), 0.0);
        if(curveArea < 2000){
            files.erase(i);
            continue;
        }
        
        //remove the temp.csv created in fileManager
        remove( (dir + "/temp.csv").c_str() );
        
        cout<<"."<<endl<<"Finding Peaks  ..";
        cout.flush();
        findPeaks(data.first,data.second, peakParams, output_dir);
        cout<<".";
        cout.flush();
        stats[count].push_back(fileManager.barcode());
        
        //need cast to cstr for remove in stdio.h
        remove( (dir + "/temp.csv").c_str() );
        
        cout.flush();
        cout<<endl<<"Deconvoluting Glow Peak  .";
        cout.flush();
        First_Order_Kinetics FOK_Model = *new First_Order_Kinetics(data,peakParams);
        stats[count].push_back(FOK_Model.glow_curve());
        stats[count].push_back(integral);
        stats[count].push_back(fileManager.temp_rate(filename));
        vector<double> peak_integral = FOK_Model.return_curve_areas();
        for(auto i = peak_integral.begin(); i != peak_integral.end();++i){
            stats[count].push_back(*i);
        }
        vector<vector<double>> returnedPeaks = FOK_Model.return_glow_curve();
        filenames.push_back(filename);
        filename = output_dir +"/"+ filename;
        fileManager.write(returnedPeaks, filename);
        cout<<"----------------------------"<<endl;
        count++;
        if(count == int(files.size()))fileManager.statistics(stats,filenames,output_dir);
    }
    return 0;
}
