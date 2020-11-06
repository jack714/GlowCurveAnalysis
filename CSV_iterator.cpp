//
//  CSV_iterator.cpp
//  GlowCurveAnalsys
//
//  Initially created by jeremy hepker on 1/13/19.
//
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
using namespace std;

//this class creates an object that takes an istream to read in csv data
template <class T>
class csv_iterator: public iterator<input_iterator_tag, T>
{
    istream * input;
    char delim;
    string value;
public:
    csv_iterator( char delm = ',' ): input( 0 ), delim( delm ) {}
    //use istream to read in data until ',' at a time
    csv_iterator( istream & in, char delm = ',' ): input( &in ), delim( delm ) { ++*this; }
    
    const T operator *() const {
        istringstream ss( value );
        T value;
        ss >> value;
        return value;
    }
    
    istream & operator ++() {
        if( !( getline( *input, value, delim ) ) )
        {
            input = 0;
        }
        return *input;
    }
    
    bool operator !=( const csv_iterator & rhs ) const {
        return input != rhs.input;
    }
};

