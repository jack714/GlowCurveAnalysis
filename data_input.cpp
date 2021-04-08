//
//  data_input.cpp
//  GlowCurveAnalsys
//
//  Created by Jack Yu UROP 2020 Fall
//

#include "data_input.hpp"

//read in user data
vector<vector<double>> input_data() {
    string input;
    vector<vector<double>> data;
    double energy_input, temp_input, count_input;
    //used to get rid of blank line
    getline(cin, input);
    while (true) {
        vector<double> peak;
        getline(cin, input);
        if (input == "done")
            break;
        //cast the string and store corresponding values
        auto index = input.find(',');
        if (index == std::string::npos) {
            cout << "Invalid input!" << endl;
            continue;
        }
        energy_input = stod(input.substr(0, index));
        auto index2 = input.find(',', index + 1);
        temp_input = stod(input.substr(index + 1, index2 - index - 1));
        count_input = stod(input.substr(index2 + 1));
        peak.push_back(energy_input);
        peak.push_back(temp_input);
        peak.push_back(count_input);
        data.push_back(peak);
    }
    return data;
}