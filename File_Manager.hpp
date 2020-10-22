//
//  File_Manager.hpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#ifndef File_Manager_hpp
#define File_Manager_hpp

#include <string>
#include <getopt.h>
#include <cmath>
#include <vector>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include <iostream>
#include <locale>
#include <cstdlib>
#include <iomanip>
#include "DataSmoothing.hpp"

class File_Manager{
private:
    std::vector<double> raw_temp_data;
    std::vector<double> raw_count_data;
    std::vector<double> heating_rate;
    int total_counts = 0;
    double maxTime = 0.0, maxTemp = 0.0, barcodeNum = 0.0;
    bool time = false;
    std::string filename,header;
public:
    
    File_Manager(std::string filename);
    File_Manager();
    
    //This function reads in the .csv file and parses the raw data into vector of coordinate pairs.
    std::pair<std::vector<double>,std::vector<double>> read();
    void statistics(std::vector<std::vector<double>> stats, std::vector<std::string> filenames, std::string dir);
    //This is a function to write the output to a new CSV file.
    void write(std::vector<std::vector<double>> glow_curves, std::string output_name);
    double temp_rate(std::string name);
    double barcode(){
        return barcodeNum;
    };
    ~File_Manager();
    
};


#endif /* File_Manager_hpp */
