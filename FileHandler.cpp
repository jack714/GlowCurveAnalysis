//	fileHandler.cpp
//	GlowCurveAnalysis
//
//	Created by Jack Yu UROP 2020 Fall

#include "FileHandler.hpp"

using namespace std;
//this needs C++17 or above
namespace fs = std::filesystem;

//open up the input dir directory and find all csv files
//return a vector of string that contains all the paths to csv files
//update output_dir to be the local directory for output folder
vector<string> handle_dir(string dir, string& output_dir) {

	int numDir = 0;
	int numCSV = 0;
	vector<string> csv;

	//check if the directory actually exist
	try {
		if (!fs::exists(dir)) {
			throw dir + " Invalid directory path!";
		}
	}
	catch (string str) {
		cerr << str << endl;
	}

	//assigned the output folder path to output_dir
    output_dir = dir + "_output";
    fs::create_directories(output_dir);

	//iteratively read in all the files/sub-directories in the "dir" directory
	for (const auto& entry : fs::recursive_directory_iterator(dir)) {
		//if find sub-directory then record the number of sub-directories found
		if (fs::is_directory(entry)) {
			numDir++;
		}
		//if it's csv/xlsx file then add to the output vector and increment number of csv found
		else {
            //cast the path and check the format
			string path = entry.path().string();
			unsigned long dot = path.find_last_of('.');
			string type = path.substr(dot + 1);
			if (type == "xlsx" || type == "csv" || type == "xls") {
				csv.push_back(path);
				numCSV++;
			}
		}
	}
    //output the details about the directory readed in
	cout << "Directories found: " << numDir << endl;
	cout << "CSV files found: " << numCSV << endl;
	cout << "Output Directory created: " << output_dir << endl;
	cout << "CSV details: " << endl;
    //for debugging purpose, shows all the files read
	//std::sort(csv.begin(), csv.end(), compare);
	for (string s : csv) {
		cout << s << endl;
	}
	return csv;
}

bool compare(string a, string b) {
	int first = stoi(a.substr(a.find_last_of("_") + 1, a.find_last_of(".")));
	int sec = stoi(b.substr(b.find_last_of("_") + 1, b.find_last_of(".")));
	return first < sec;
}