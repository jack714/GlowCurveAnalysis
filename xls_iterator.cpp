//
//  xls_iterator.cpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu UROP 2020 Fall
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
using namespace std;

//this class creates an object that takes an istream to read in xls data
template <class T>
class xls_iterator : public iterator<input_iterator_tag, T>
{
    istream* input;
    char delim;
    string value;
public:
    xls_iterator(char delm = '\t') : input(0), delim(delm) {}
    //use istream to read in data one at a time
    xls_iterator(istream& in, char delm = '\t') : input(&in), delim(delm) { ++* this; }

    const T operator *() const {
        istringstream ss(value);
        T value;
        ss >> value;
        return value;
    }

    istream& operator ++() {
        if (!(getline(*input, value, delim)))
        {
            input = 0;
        }
        return *input;
    }

    bool operator !=(const xls_iterator& rhs) const {
        return input != rhs.input;
    }
};