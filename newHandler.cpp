#include "newHandler.hpp"

using namespace std;
//this needs C++17 or above
namespace fs = std::filesystem;

//open up the input dir directory and find all csv files
vector<string> handle_dir(string dir) {
	//this is the path for my machine
	//string pathing = "C:/Users/Wenji/Desktop/GlowCurveAnalsys-master_new/GlowCurveAnalsys-master_new/";
    //string pathing = "/Users/rhellab/Desktop/GlowCurveAnalysis-master_new/GlowCurveAnalsys-master_new/";
	//dir = pathing + dir;

	int numDir = 0;
	int numCSV = 0;
	vector<string> csv;

	//check if the directory actually exist, need full path for this
	//improvement can be finding file in a set directory aka GCA folder
	try {
		if (!fs::exists(dir)) {
			throw dir + " Invalid directory path!";
		}
	}
	catch (string str) {
		cerr << str << endl;
		//return -1;
	}

	//create the output folder, in Visual Studio the folder is created at the same folder with the code
	//improvement can be setting the new directory to be placed in the main folder
	unsigned long last = dir.find_last_of('/');
	string folder = dir.substr(last + 1);
	string output = folder + "_output";
	fs::create_directories(output);

	//iteratively read in all the files/sub-directories in the "dir" directory
	for (const auto& entry : fs::recursive_directory_iterator(dir)) {
		//if find sub-directory then record it
		if (fs::is_directory(entry)) {
			numDir++;
		}
		//if it's csv file then add to the output vector and increment number of csv found
		else {
			string path = entry.path().string();
			unsigned long dot = path.find_last_of('.');
			string type = path.substr(dot + 1);
			if (type == "xlsx" || type == "csv") {
				csv.push_back(path);
				numCSV++;
			}
		}
	}
	cout << "Directories found: " << numDir << endl;
	cout << "CSV files found: " << numCSV << endl;
	cout << "Output Directory created: " << output << endl;
	cout << "CSV details: " << endl;
	for (string s : csv) {
		cout << s << endl;
	}
	return csv;
}
