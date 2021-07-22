#include <iostream>
#include <getopt.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <set>
#include "FileHandler.hpp"
using namespace std;
namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    srand(100);
    string dir;
    string output_dir = "";
    cout << "Please enter the full path to directory containing csv formatted emission spectra:" << endl;
    cin >> dir;
    vector<string> files;
    files = handle_dir(dir, output_dir);
    string test = output_dir + "_test";
    string train = output_dir + "train";
    string verify = output_dir + "verify";
    fs::create_directories(test);
    fs::create_directories(train);
    fs::create_directories(verify);
    //choose train set, size is 80
    set<int> s;
    int size = files.size();
    for (int i = 0; i < 80; i++) {
        int r = rand() % size;
        if (s.find(r) != s.end())
            i--;
        else {
            fs::copy(files[r], fs::path(train));
            s.insert(r);
        }
    }
    //choose test set
    int test_size = size * 2 / 5;
    for (int i = 0; i < test_size; i++) {
        int r = rand() % size;
        if (s.find(r) != s.end())
            i--;
        else {
            fs::copy(files[r], fs::path(test));
            s.insert(r);
        }
    }
    //choose verify
    for (int i = 0; i < size - 80 - test_size; i++) {
        int r = rand() % size;
        if (s.find(r) != s.end())
            i--;
        else {
            fs::copy(files[r], fs::path(verify));
            s.insert(r);
        }
    }
}