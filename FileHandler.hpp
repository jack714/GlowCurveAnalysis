//	fileHandler.hpp
//	GlowCurveAnalysis
//
//	Created by Jack Yu UROP 2020 Fall

#ifndef newHandler_hpp
#define newHandler_hpp
#include <iostream>
#include <string>
#include <vector>
//filesystem is available in C++17, please specify this c++ version
//if a compiler like visual studio is used. If using makefile then no changes needed
#include <filesystem>
using namespace std;

//open up the dir directory and store all the csv file names
//update output_dir to be the local directory for output folder
std::vector<std::string> handle_dir(std::string dir, std::string& output_dir);
//sort the files in order
bool compare(std::string a, std::string b);
#endif /* newHandler_hpp */


