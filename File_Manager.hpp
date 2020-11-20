//
//  File_Manager.hpp
//  GlowCurveAnalsys
//
//  Initially created by jeremy hepker on 1/27/19.
//
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//

#ifndef File_Manager_hpp
#define File_Manager_hpp

#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include <stdio.h>
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
    std::string filename, header, output_dir;
public:
    
    File_Manager(std::string filename, std::string out);
    File_Manager();
    
    //This function reads in the .csv file and parses the raw data into vector of coordinate pairs.
    std::pair<std::vector<double>,std::vector<double>> read(double& time);
    void statistics(std::vector<std::vector<double>> stats, std::vector<std::string> filenames, std::string dir, int count);
    //This is a function to write the output to a new CSV file.
    void write(std::vector<std::vector<double>> glow_curves, std::string output_name);
    double temp_rate(std::string name);
    double barcode(){
        return barcodeNum;
    };
    ~File_Manager();
    
};


#endif /* File_Manager_hpp */
