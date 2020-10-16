//
//  File_Manager.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "File_Manager.hpp"
#include "CSV_iterator.cpp"

File_Manager::File_Manager(std::string given_filename):filename(given_filename){};

//This function reads in the .csv file and parses the data into std::vector of coordinate pairs.
pair<std::vector<double>,std::vector<double>> File_Manager::read(){
    //Open and test the user input file.
    std::string line, catagories;
    vector<double> tempData, countData;
    size_t path = filename.find_last_of("/\\");
    std::string temp_path = filename.substr(0,path+1);
    temp_path += "temp.csv";
    //check if filename is a valid path
    ifstream file(filename);
    if(!file.is_open()){
        cerr<<"Error opening file: "<<filename<<endl;
        exit(1);
    }
    //get rid of the title content
    while(true){
        if(line.find("Time (") != std::string::npos) time =true;
        if(line.find("Count")==std::string::npos){
            getline(file, line,'\n');
        }else{
            getline(file, line,'\n');
            break;
        }
    }
    ofstream temp_file;
    temp_file.open(temp_path);
    if(!temp_file.is_open()){
        cerr<<"Error opening file: "<<filename<<endl;
        exit(1);
    }
    std::vector<std::string> temps;
    while(getline(file, line,'\n')){
        if(line.find('\r') != std::string::npos) line.pop_back();
        if(line.back() != ',') temps.push_back(line + ",\n");
        else temps.push_back(line+ "\n");
    }
    temps.back().pop_back();
    temps.back().pop_back();
    for(auto j = temps.begin(); j != temps.end(); ++j){
        temp_file<<*j;
    }
    temp_file.close();
    file.close();
    file.open(temp_path);
    getline(file, line);
    getline(file, line);
    std::stringstream ss;
    ss << line;
    int count = 0;
    while(getline(ss, line, ',')){
        count++;
        if(line.find('N') != std::string::npos){
            line.erase(0,1);
            barcodeNum = stoi(line);
        }
    }
    auto i = csv_iterator<std::string>( file );
    bool two = false;
    while(file){
        for(int j=0;j<(count-2);++j) ++i;
        raw_temp_data.push_back(stod(*i));
        ++i;
        raw_count_data.push_back(stod(*i));
        if(!two && raw_count_data.back() >2)two = true;
        ++i;
        if(file.eof()) break;
        if(!two){
            raw_temp_data.pop_back();
            raw_count_data.pop_back();
        }
    }
    file.close();
    if((raw_temp_data.size()%2) == 0){
        raw_temp_data.pop_back();
        raw_count_data.pop_back();
    }
    
    // *************** //
    // *************** //
    // *************** //
    /*
    // average first 500 points
    double total = 0;
    int num = int( raw_temp_data.size( ) );
    for ( int i = 0 ; i < num ; i++ )
    {
        total += raw_temp_data[i];
    }
    double average = total / num;
    */
    // *************** //
    // *************** //
    // *************** //
    
    int n1 = 0, n2 = 100;
    double x1 = 0.0, x2 = 0.0, y1 = 0.0, y2 = 0.0;
    x1 = raw_count_data[n1];
    y1 = raw_temp_data[n1];
    x2 = raw_count_data[n2];
    y2 = raw_temp_data[n2];
    double m = ( y2 - y1 ) / ( x2 - x1 ), b = y1 - ( m * x1 );
    double linear_approximation = 0.0;
    
    for ( int i = 0 ; i < int( raw_temp_data.size( ) ) ; i++ )
    {
        linear_approximation = m * raw_count_data[i] + b;
        if ( raw_temp_data[i] < linear_approximation )
        {
            raw_temp_data[i] = 0;
        }
        else
        {
            raw_temp_data[i] -= linear_approximation;
        }
    }
    
    for(int i = 0; i < 5; ++i) dataSmooth(raw_temp_data, raw_count_data);
    
    // *************** //
    // *************** //
    // *************** //
    /*
    for ( int i = 0 ; i < num ; i++ )
    {
        if ( raw_temp_data[i] < average )
        {
            raw_temp_data[i] = 0;
        }
        else
        {
            raw_temp_data[i] -= average;
        }
     }
     */
    // *************** //
    // *************** //
    // *************** //
    
    return make_pair(raw_temp_data, raw_count_data);
}

//This is a function to write the output to a new CSV file.
void File_Manager::write(std::vector<std::vector<double>> glow_curves, std::string output_name){
    ofstream file;
    output_name += "_output.csv";
    file.open(output_name);
    if(!file.is_open()){
        cerr<<"Could not open output file : "<<output_name<<endl;
        exit(1);
    }
    file<<"temp,count_original";
    for(int j = 0; j<int(glow_curves.size());++j){
        std::string ster = "count_" + std::to_string(j);
        file<<","<<ster;
    }
    file<<",\n";
    file.setf(ios_base::fixed);
    file<<setprecision(5);
    for(int i = 0; i<int(raw_temp_data.size());++i){
        file << raw_temp_data[i]<<",";
        file << raw_count_data[i];
        for(int j = 0; j<int(glow_curves.size());++j){
            file<<","<<double(glow_curves[j][i]);
        }
        file<<",\n";
    }
    cout<<"Output File : "<<output_name<<endl;
    file.close();
}

double File_Manager::temp_rate(std::string name){
    if(isdigit(name[0])){
        int j = 0;
        while(isdigit(name[j])){
            ++j;
            if(name[j] == '.') ++j;
        }
        std::string temp = name.substr(0,j);
        return(stod(temp));
    }
    auto max = max_element(raw_temp_data.begin(), raw_temp_data.end());
    int pos = int(max - raw_temp_data.begin());
    return (raw_temp_data[pos]-raw_temp_data[8])/(int(max - raw_temp_data.begin())/2);
}

void File_Manager::statistics(std::vector<std::vector<double>> stats, std::vector<std::string> filenames, std::string dir){
    std::ofstream myfile;
    dir += "/statistics.csv";
    myfile.open(dir);
    if(!myfile.is_open()){
        std::cerr<<"Could not open output file : statistics.csv"<<std::endl;
        exit(1);
    }
    myfile<<"Filename, Barcode, Figure of Merit, Total Curve Area, Heating Rate (C/s), Peak Area(s)\n";
    int count = 0; 
    for(auto i = stats.begin(); i != stats.end();++i){
        myfile<<filenames[count++]<<",";
        for(auto j = i->begin(); j != i->end();++j){
            myfile<<*j<<",";
        }
        myfile<<"\n";
    }
    count = 0;
    myfile.close();
}

File_Manager::~File_Manager()
{
    raw_temp_data.clear();
    raw_count_data.clear();
    heating_rate.clear();
}
