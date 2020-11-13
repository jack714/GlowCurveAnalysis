//
//  File_Manager.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//

#include "File_Manager.hpp"
#include "CSV_iterator.cpp"
#include "background_subtraction.hpp"

File_Manager::File_Manager(std::string given_filename):filename(given_filename){};

//This function reads in the .csv file and parses the raw data into std::vector of coordinate pairs.
pair<std::vector<double>,std::vector<double>> File_Manager::read(){
    //Open and test the user input file.
    std::string line;
    vector<double> tempData, countData;
    //create a temporary csv file to store data
    size_t path = filename.find_last_of("/\\");
    std::string temp_path = filename.substr(0,path+1);
    temp_path += "temp.csv";
    //open and check if filename is a valid path
    ifstream file(filename);
    if(!file.is_open()){
        cerr<<"Error opening file: "<<filename<<endl;
        exit(1);
    }
    //get rid of the title content, this will also get rid of the first line of the actual data
    while(true){
        if(line.find("Time (") != std::string::npos) time =true;
        if(line.find("Count") == std::string::npos){
            getline(file, line,'\n');
        }else{
            getline(file, line,'\n');
            break;
        }
    }
    //open the temp.csv file
    ofstream temp_file;
    temp_file.open(temp_path);
    if(!temp_file.is_open()){
        cerr<< "Error opening file: " << filename << endl;
        exit(1);
    }
    //reading in each line of data to vector temps
    std::vector<std::string> temps;
    while(getline(file, line,'\n')){
        if(line.find('\r') != std::string::npos)
            line.pop_back();
        if(line.back() != ',')
            temps.push_back(line + ",\n");
        else
            temps.push_back(line + "\n");
    }
    temps.back().pop_back();
    temps.back().pop_back();
    //copy each line in temp to the temp.csv file
    for(auto j = temps.begin(); j != temps.end(); ++j){
        temp_file << *j;
    }
    temp_file.close();
    file.close();
    
    //open the temporary csv file and process data
    file.open(temp_path);
    getline(file, line);
    getline(file, line);
    std::stringstream ss;
    ss << line;
    int count = 0;
    //process the first line of data, for each time reading in until ','
    //use count to record the number of ',' until temperature data
    while(getline(ss, line, ',')){
        count++;
        //if there's a column for barcode, this would get executed and barcodeNum would record the barcode
        if(line.find('N') != std::string::npos){
            line.erase(0,1);
            barcodeNum = stoi(line);
        }
    }
    //construct an csv_iterator object from csv_iterator.cpp to read in data from temp.csv
    auto i = csv_iterator<std::string>( file );
    /*
     I feel this can be improved by getting rid of two
     */
    //bool two = false;
    while(file){
        //skip the first count-2 sets of data separtaed by ',' so the next data read in is temperature
        for(int j = 0; j < (count - 2); ++j)
            ++i;
        //push the temperature data to raw_temp_data
        raw_temp_data.push_back(stod(*i));
        //read in data until next ','
        ++i;
        //push count data to raw_count_data vector
        raw_count_data.push_back(stod(*i));
        //if(raw_count_data.back() <= 2) {
        //    raw_temp_data.pop_back();
        //    raw_count_data.pop_back();
        //}
        //if(!two && raw_count_data.back() > 2)
        //    two = true;
        //go to next data by increment i again
        ++i;
        if(file.eof())
            break;
        //get rid of the data with count smaller or equal to 2
        //if(!two){
        //    raw_temp_data.pop_back();
        //    raw_count_data.pop_back();
        //}
    }
    // BACKGROUND_SUBTRACTION
    bg_subtract(raw_temp_data, raw_count_data);
    //get rid of count < 2 data
    int index = 0;
    while (raw_count_data[index] < 2) {
        if (!raw_count_data.empty()) {
            raw_count_data.erase(raw_count_data.begin() + index);
            raw_temp_data.erase(raw_temp_data.begin() + index);
        }
        if (raw_count_data.empty()) {
            break;
        }
        //raw_temp_data.erase(raw_temp_data.begin() + index);
    }
    file.close();
    //if the size of the two vector is even then remove the last data to make the size odd
    if((raw_temp_data.size() % 2) == 0 && !raw_temp_data.empty()){
        raw_temp_data.pop_back();
        raw_count_data.pop_back();
    }
    
    return make_pair(raw_temp_data, raw_count_data);
}

//This is a function to write the output to a new CSV file.
void File_Manager::write(std::vector<std::vector<double>> glow_curves, std::string output_name){
    ofstream file;
    output_name += "_output.csv";
    file.open(output_name);
    if(!file.is_open()){
        cerr << "Could not open output file : "<<output_name<<endl;
        exit(1);
    }
    //output the titles including temp, original count, and new counts for peak fitting
    file << "temp,count_original";
    for(int j = 0; j < int(glow_curves.size()); ++j){
        std::string ster = "count_" + std::to_string(j);
        file << ","<<ster;
    }
    file << ",\n";
    file.setf(ios_base::fixed);
    file << setprecision(5);
    //output every temperature, original count, and FOK count value under each peak curve fit
    //for(int i = 0; i < int(raw_temp_data.size()); ++i){
    for (int i = 0; i < int(glow_curves[0].size()); ++i) {
        file << raw_temp_data[i] << ",";
        file << raw_count_data[i];
        for(int j = 0; j<int(glow_curves.size());++j){
            file << "," << double(glow_curves[j][i]);
        }
        file << ",\n";
    }
    //output the file name to make sure the address is correct
    cout << "Output File : " << output_name << endl;
    file.close();
}

double File_Manager::temp_rate(std::string name){
    //cast the name and assign the numeric portion which is temperature rate to temp
    if(isdigit(name[0])){
        int j = 0;
        while(isdigit(name[j])){
            ++j;
            if(name[j] == '.') ++j;
        }
        std::string temp = name.substr(0,j);
        return(stod(temp));
    }
    //a way to compute a temperature rate if the file name doesn't contain the rate
    auto max = max_element(raw_temp_data.begin(), raw_temp_data.end());
    int pos = int(max - raw_temp_data.begin());
    return (raw_temp_data[pos] - raw_temp_data[8]) / (int(max - raw_temp_data.begin()) / 2);
}

//output an overview of all processed fittings 
void File_Manager::statistics(std::vector<std::vector<double>> stats, std::vector<std::string> filenames, std::string dir, int num){
    std::ofstream myfile;
    dir += "/statistics.csv";
    myfile.open(dir);
    if(!myfile.is_open()){
        std::cerr<<"Could not open output file : statistics.csv"<<std::endl;
        exit(1);
    }
    myfile<<"Filename, Barcode, Figure of Merit, Total Curve Area, Heating Rate (C/s), Peak Area(s)\n";
    int count = 0; 
    for(int i = 0; i < num;++i){
        myfile << filenames[count++] << ",";
        for(auto j = stats[i].begin(); j != stats[i].end(); ++j){
            myfile << *j << ",";
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
