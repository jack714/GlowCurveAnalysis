//
//  CSV_iterator.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/13/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <sstream>
using namespace std;

template <class T>
class csv_iterator: public iterator<input_iterator_tag, T>
{
    istream * input;
    char delim;
    string value;
public:
    csv_iterator( char delm = ',' ): input( 0 ), delim( delm ) {}
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

