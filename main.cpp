//
//  main.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include <iostream>
#include <getopt.h>
#include <string>
#include <locale>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include "File_Manager.hpp"
//#include "batch_handler.hpp"
#include "newHandler.hpp"
#include "smartPeakDetect.hpp"
#include "DataSmoothing.hpp"
#include "Levenberg–Marquardt.hpp"
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;

int main() {
    string dir,start ="n";
    vector<string> filenames,files;
    while(start =="n" || start =="N"){
        cout<<"Please enter the full path to directory containing csv formatted emission spectra:"<<endl;
        cin>>dir;
        //string dir = argv[1];
        if(dir.back() == '/') 
            dir.pop_back();
        //files = batch_handler(dir);
        files = handle_dir(dir);
        cout<<"Is this correct and would you like to start processing (y/n)?"<<endl;
        cin>>start;
    }
    vector<vector<double>> stats(files.size(), vector<double>(0,0.0));
    int count = 0;
    auto i = files.begin();
    string output_dir = *i++;
    for(; i != files.end(); ++i){
        vector<vector<double>> peakParams;
        double integral= 0;
        if(i->find("temp.csv") !=string::npos){
            files.erase(i);
            continue;
        }
        cout<<"----------------------------"<<endl<<"Processing: ";
        string filename = i->substr((i->find_last_of("/\\"))+1);
        cout<<filename<<" ("<<count+1<<" of "<<files.size()-1<<")"<<endl<<"Reading in File  .";
        cout.flush();
        File_Manager fileManager = *new File_Manager(*i);
        cout<<".";
        cout.flush();
        pair<vector<double>, vector<double>> data = fileManager.read();
        const double curveArea = accumulate(data.second.begin(), data.second.end(), 0.0);
        if(curveArea < 2000){
            files.erase(i);
            continue;
        }
        //remove( dir + "/temp.csv" );
        cout<<"."<<endl<<"Finding Peaks  ..";
        cout.flush();
        findPeaks(data.first,data.second, peakParams);
        cout<<".";
        cout.flush();
        stats[count].push_back(fileManager.barcode());
        //remove( dir + "/temp.csv" );
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
