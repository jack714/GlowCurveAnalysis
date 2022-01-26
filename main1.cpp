//
//  main.cpp
//  GlowCurveAnalsys
//
//  Initially created by jeremy hepker on 1/9/19.
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//

#include <iostream>
#include <getopt.h>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include "File_Manager.hpp"
#include "FileHandler.hpp"
#include "smartPeakDetect.hpp"
#include "DataSmoothing.hpp"
#include "Levenberg_Marquardt.hpp"
#include "quick_half_max.hpp"
#include "data_input.hpp"
#include "remove_spike.hpp"
#include "Savitzky_Golay.hpp"
#include "background_remove.hpp"
#include <chrono>
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
//if unistd.h can't be included then comment out the else clause
//#else
//#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;

double func(const double input, const vector<double> params) {
    double T = 0.0;
    double I_t = 0.0;
    double k = .000086173303;
    double energy = params[0];
    double Tm = params[1] + 273.15;
    double dm = (2.0 * k * (Tm)) / energy;
    double Im = params[2];
    T = double(input + 273.15);
    I_t = Im * exp(1.0 + (energy / (k * T)) * ((T - Tm) / Tm) - ((T * T) / (Tm * Tm)) * exp((energy / (k * T)) * ((T - Tm) / Tm)) * (1.0 - ((2.0 * k * T) / energy)) - dm);
    if (isinf(I_t)) {
        cout << T << " " << Tm << " " << Im << " " << dm << endl;
    }
    return I_t;
}

double quick_fom(vector<double>& temp, vector<double>& count, vector<vector<double>>& peak_param) {
    vector<double> total_curve(temp.size());
    double orig_integral = 0.0;
    vector<double> total_fit(count.size());
    for (int d = 0; d < int(count.size()); d++) {
        double fit = 0.0;
        for (int e = 0; e < int(peak_param.size()); e++) {
            fit += func(temp[d], peak_param[e]);
        }
        orig_integral += fit;
        total_fit[d] = fit;
    }
    double fom = 0.0;
    for (int f = 0; f < int(count.size()); ++f) {
        fom += abs(count[f] - total_fit[f]) / orig_integral;
    }
    return fom;
}

void output_details(string output_dir, string filename, vector<double>& temp, vector<double> count, vector<vector<double>>& peak_300, vector<vector<double>>& orig_peak_300,
    vector<vector<double>>& peak_400, vector<vector<double>>& orig_peak_400, vector<vector<double>>& peak_900, vector<vector<double>>& orig_peak_900,
    vector<vector<double>>& peak_100, vector<vector<double>>& orig_peak_100, vector<vector<double>>& peak_200, vector<vector<double>>& orig_peak_200) {
    ofstream file2;
    string path = output_dir + "/300_" + filename;
    file2.open(path);
    vector<vector<double>> gd_300(peak_300.size(), vector<double>(count.size(), 0));
    vector<vector<double>> orig_300(peak_300.size(), vector<double>(count.size(), 0));
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(orig_peak_300.size()); ++x) {
            double out = quickFok(count[i], orig_peak_300[x]);
            orig_300[x][i] = out;
        }
    }
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(peak_300.size()); ++x) {
            double out = quickFok(temp[i], peak_300[x]);
            gd_300[x][i] = out;
        }
    }
    file2 << "temp, after removal, first, sec, third, forth, fifth, sixth, seventh, eighth, GDfirst, GDsec, GDthird, GDforth, GDfifth, GDsixth, GDseventh, GDeighth";
    file2 << ",\n";
    for (int i = 0; i < int(temp.size()); i++) {
        file2 << temp[i] << ",";
        file2 << count[i] << ",";
        for (int j = 0; j < int(orig_peak_300.size()); j++) {
            file2 << orig_300[j][i] << ",";
        }
        for (int j = 0; j < int(peak_300.size()) - 1; j++) {
            file2 << gd_300[j][i] << ",";
        }
        file2 << gd_300[peak_300.size() - 1][i];
        file2 << ",\n";

    }
    file2.close();

    ofstream file3;
    path = output_dir + "/400_" + filename;
    file3.open(path);
    vector<vector<double>> gd_400(peak_400.size(), vector<double>(count.size(), 0));
    vector<vector<double>> orig_400(peak_400.size(), vector<double>(count.size(), 0));
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(orig_peak_400.size()); ++x) {
            double out = quickFok(temp[i], orig_peak_400[x]);
            orig_400[x][i] = out;
        }
    }
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(peak_400.size()); ++x) {
            double out = quickFok(temp[i], peak_400[x]);
            gd_400[x][i] = out;
        }
    }
    file3 << "temp, after removal, first, sec, third, GDfirst, GDsec, GDthird";
    file3 << ",\n";
    for (int i = 0; i < int(temp.size()); i++) {
        file3 << temp[i] << ",";
        file3 << count[i] << ",";
        for (int j = 0; j < int(orig_peak_400.size()); j++) {
            file3 << orig_400[j][i] << ",";
        }
        for (int j = 0; j < int(peak_400.size()) - 1; j++) {
            file3 << gd_400[j][i] << ",";
        }
        file3 << gd_400[peak_400.size() - 1][i];
        file3 << ",\n";

    }
    file3.close();

    ofstream file4;
    path = output_dir + "/900_" + filename;
    file4.open(path);
    vector<vector<double>> gd_900(peak_900.size(), vector<double>(count.size(), 0));
    vector<vector<double>> orig_900(peak_900.size(), vector<double>(count.size(), 0));
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(orig_peak_900.size()); ++x) {
            double out = quickFok(temp[i], orig_peak_900[x]);
            orig_900[x][i] = out;
        }
    }
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(peak_900.size()); ++x) {
            double out = quickFok(temp[i], peak_900[x]);
            gd_900[x][i] = out;
        }
    }
    file4 << "temp, after removal, first, sec, third, fourth, fifth, sixth, seventh, GDfirst, GDsec, GDthird, GDfourth, GDfifth, GDsixth, GDseventh";
    file4 << ",\n";
    for (int i = 0; i < int(temp.size()); i++) {
        file4 << temp[i] << ",";
        file4 << count[i] << ",";
        for (int j = 0; j < int(orig_peak_900.size()); j++) {
            file4 << orig_900[j][i] << ",";
        }
        for (int j = 0; j < int(peak_900.size()) - 1; j++) {
            file4 << gd_900[j][i] << ",";
        }
        file4 << gd_900[peak_900.size() - 1][i];
        file4 << ",\n";

    }
    file4.close();
    ofstream file5;
    path = output_dir + "/100_" + filename;
    file5.open(path);
    vector<vector<double>> gd_100(peak_100.size(), vector<double>(count.size(), 0));
    vector<vector<double>> orig_100(peak_100.size(), vector<double>(count.size(), 0));
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(orig_peak_100.size()); ++x) {
            double out = quickFok(temp[i], orig_peak_100[x]);
            orig_100[x][i] = out;
        }
    }
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(peak_100.size()); ++x) {
            double out = quickFok(temp[i], peak_100[x]);
            gd_100[x][i] = out;
        }
    }
    file5 << "temp, after removal, first, sec, third, forth, GDfirst, GDsec, GDthird, GDforth";
    file5 << ",\n";
    for (int i = 0; i < int(temp.size()); i++) {
        file5 << temp[i] << ",";
        file5 << count[i] << ",";
        for (int j = 0; j < int(orig_peak_100.size()); j++) {
            file5 << orig_100[j][i] << ",";
        }
        for (int j = 0; j < int(peak_100.size()) - 1; j++) {
            file5 << gd_100[j][i] << ",";
        }
        file5 << gd_100[peak_100.size() - 1][i];
        file5 << ",\n";

    }
    file5.close();
    ofstream file6;
    path = output_dir + "/200_" + filename;
    file6.open(path);
    vector<vector<double>> gd_200(peak_200.size(), vector<double>(count.size(), 0));
    vector<vector<double>> orig_200(peak_200.size(), vector<double>(count.size(), 0));
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(orig_peak_200.size()); ++x) {
            double out = quickFok(temp[i], orig_peak_200[x]);
            orig_200[x][i] = out;
        }
    }
    for (int i = 0; i < int(temp.size()); ++i) {
        double output = 0.0;
        for (int x = 0; x < int(peak_200.size()); ++x) {
            double out = quickFok(temp[i], peak_200[x]);
            gd_200[x][i] = out;
        }
    }
    file6 << "temp, after removal, first, sec, third, forth, fifth, sixth, seventh, eighth, nineth, GDfirst, GDsec, GDthird, GDforth, GDfifth, GDsixth, GDseventh, GDeighth, GDninth";
    file6 << ",\n";
    for (int i = 0; i < int(temp.size()); i++) {
        file6 << temp[i] << ",";
        file6 << count[i] << ",";
        for (int j = 0; j < int(orig_peak_200.size()); j++) {
            file6 << orig_200[j][i] << ",";
        }
        for (int j = 0; j < int(peak_200.size()) - 1; j++) {
            file6 << gd_200[j][i] << ",";
        }
        file6 << gd_200[peak_200.size() - 1][i];
        file6 << ",\n";

    }
    file6.close();
}

void output_counts(vector<double>& orig_counts, vector<double>& smoothed_count, vector<double>& temp, string filename) {
    ofstream file1;
    string path = "C:/Users/jack0/Desktop/" + filename + ".csv";
    file1.open(path);
    file1 << "temp, orig_count, current_count," << endl;
    for (int i = 0; i < int(orig_counts.size()); i++) {
        file1 << temp[i] << ", " << orig_counts[i] << ", " << smoothed_count[i] << "," << endl;
    }
}

void gd(const vector<double>& temp, const vector<double>& curve, vector<vector<double>>& peakParams, double& FOM, ofstream& f1, bool& one) {
    //cout << peakParams.size() << endl;
    //temperary vector to store peak data
    vector<vector<double>> temp_params = peakParams;
    //temperary vector to store accumulated fitted count
    vector<double> temp_output(curve.size(), 0.0);
    int curveSize = int(curve.size());
    int peakNum = int(peakParams.size());
    double k = .000086173303;
    double rate1 = 0.0000001;
    //double rate1 = 0.0000000047;
    double rate2 = 0.00002;
    //double rate2 = 0.000001;
    //double rate3 = 0.001;
    double rate3 = 0.00001;
    //double rate3 = 0.0000065;
    double current_FOM = FOM;
    int iteration = 0;
    bool check = true;
    while (iteration < 500 && check) {
        vector<vector<double>> update(peakNum, vector<double>(3, 0.0));
        vector<vector<double>> curves(temp_params.size(), vector<double>(curve.size()));
        vector<double> totalCurve(curve.size());
        for (int i = 0; i < int(curve.size()); ++i) {
            double out;
            double accum = 0.0;
            for (int x = 0; x < int(temp_params.size()); ++x) {
                out = quickFok(temp[i], temp_params[x]);
                curves[x][i] = out;
                accum += out;
            }
            totalCurve[i] = accum;
        }
        //use FWHM to find the left and right half max points
        for (int b = 0; b < peakNum; b++) {
            int TL = temp_params[b][3];
            int TR = temp_params[b][5];
            double energy = temp_params[b][0];
            double Tm = temp_params[b][1];
            double Im = temp_params[b][2];
            for (int index = TL; index < TR + 1; index++) {
            //for (int index = 0; index < int(curve.size()); index++) {
                double y = curve[index] - totalCurve[index] + curves[b][index];
                if (y < 0)
                    y = 0;
                double T = temp[index];
                double deriv_E = -2.0 * Im * (-(2.0 * k * T * T * T * exp((energy * (T - Tm)) / (Tm * k * T))) / (energy * energy * Tm * Tm) + (2.0 * Tm * k) /
                    (energy * energy) - (T * (T - Tm) * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm * Tm * k) + (T - Tm) /
                    (Tm * k * T)) * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) /
                    (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) /
                        (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0));
                double deriv_Tm = -2.0 * Im * ((2.0 * T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm * Tm) - (T * T * (1.0 - (2.0 * k * T) / energy)
                    * exp((energy * (T - Tm)) / (Tm * k * T)) * (-(energy * (T - Tm)) / (Tm * Tm * k * T) - energy / (Tm * k * T))) / (Tm * Tm) - (energy * (T - Tm)) / (Tm * Tm * k * T) -
                    energy / (Tm * k * T) - (2.0 * k) / energy) * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) -
                    (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) -
                        (2.0 * Tm * k) / energy + 1.0));
                double deriv_Im = -2.0 * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) /
                    (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) /
                    (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0));
                update[b][0] += rate1 * deriv_E;
                update[b][1] += rate2 * deriv_Tm;
                update[b][2] += rate3 * deriv_Im;
            }
        }
        //vector<double> change(peakParams.size(), -1000);
        //for (int i = 0; i < peakNum; i++) {
        //    cout << update[i][0] << "   " << update[i][1] << "   " << update[i][2] << endl;
        //    if (abs(update[i][0]) > change[0])
        //        change[0] = abs(update[i][0]);
        //    if (abs(update[i][1]) > change[1])
        //        change[1] = abs(update[i][1]);
        //    if (abs(update[i][2]) > change[2])
        //        change[2] = abs(update[i][2]);
        //    //if (abs(update[i][2]) > change[2])
        //    //    change[2] = update[i][2];
        //    if (abs(change[0]) < 0.00001 || abs(change[1]) < 0.0001 || abs(change[2]) < 0.0001) {
        //    //if (abs(change[0]) < 0.00001 || change[2] > 0) {
        //        check = false;
        //    }
        //}
        //cout << endl;
        //apply the change to the peak paramter
        for (int c = 0; c < peakNum; c++) {
            temp_params[c][0] -= update[c][0];
            temp_params[c][1] -= update[c][1];
            temp_params[c][2] -= update[c][2];
        }
        //for (int j = 0; j < peakNum; j++) {
        //    cout << update[j][0] << ", " << update[j][1] << ", " << update[j][2] << endl;
        //    cout << temp_params[j][0] << ", " << temp_params[j][1] << ", " << temp_params[j][2] << endl;
        //}
        double new_integral = 0.0;
        vector<double> new_output(curve.size(), 0.0);
        for (int d = 0; d < curveSize; d++) {
            double new_fit = 0.0;
            for (int e = 0; e < peakNum; e++) {
                new_fit += func(temp[d], temp_params[e]);
            }
            new_integral += new_fit;
            new_output[d] = new_fit;
        }
        double new_fom = 0.0;
        for (int f = 0; f < curveSize; ++f) {
            new_fom += abs(curve[f] - new_output[f]) / new_integral;
        }
        if (new_fom > current_FOM)
            check = false;
        else
            current_FOM = new_fom;
        cout << ".";
        cout.flush();
        iteration++;
    }
    if (iteration == 1)
        one = true;
    FOM = current_FOM;
    peakParams = temp_params;
    //cout << "GD ran " << iteration << " times" << endl;
    //cout << "the new FOM is: " << FOM;
    f1 << " new fom: " << FOM << " " << iteration << " times" << endl;
}

void remove_tail(vector<double>& temp, vector<double>& curve, double intensity) {
    //find index of max point
    int max_index = -1;
    for (int i = 0; i < int(curve.size()); i++) {
        if (curve[i] == intensity)
            max_index = i;
    }
    vector<double> firstDir(temp.size(), 0.0);
    //calculate first derivative of count data using function from SmartPeakDetect
    firstDeriv(temp, curve, firstDir);
    int cut = -1;
    for (int i = max_index + 5; i < int(temp.size()); i++) {
        if (firstDir[i] > 0) {
            cut = i;
            break;
        }
    }
    cout << cut;
    temp.resize(cut);
    curve.resize(cut);
}

bool out_bounds(vector<vector<double>>& peakParams, vector<double> temp) {
    for (int i = 0; i < peakParams.size(); i++) {
        if (peakParams[i][1] > temp[temp.size() - 1])
            return true;
    }
    return false;
}

int gd_types(const vector<double>& temp, const vector<double>& curve, vector<vector<double>>& peakParams, double& FOM, int type, double max_intensity, double& original_fom, vector<vector<double>>& orig_peak) {
    //find index of max point
    int max_index = -1;
    for (int i = 0; i < int(curve.size()); i++) {
        if (curve[i] == max_intensity)
            max_index = i;
    }
    //cout << temp[max_index] << endl;
    vector<vector<double>> change_range;
    double intensity_coeff;
    double index_coff;
    double change_range_coeff;
    double best_fom = -1;
    if (type == 100) {
        change_range_coeff = 0.4;
        change_range = { {0.1, 0.05, change_range_coeff}, {0.15, 0.04, change_range_coeff}, {0.16, 0.03, change_range_coeff}, {0.15, 0.06, change_range_coeff * 1.5} };
        index_coff = temp[max_index] - 237;
        if (abs(index_coff) > 30)
            index_coff = copysign(30.0, index_coff);

        peakParams = { { 1.56, temp[max_index] - 94, max_intensity / 5, 0, 0, 0 },
                    { 1.67, temp[max_index] - 54, max_intensity / 3, 0, 0, 0},
                    { 1.69, temp[max_index] - 26, max_intensity / 2, 0, 0, 0},
                    { 2.04, temp[max_index], max_intensity, 0, 0, 0} };
        if (out_bounds(peakParams, temp)) {
            FOM = 1;
            return 1000;
        }
    }
    //else if (type == 101) {
    //    change_range_coeff = 0.2;
    //    change_range = { {0.1, 0.05, change_range_coeff}, {0.1, 0.05, change_range_coeff}, {0.1, 0.05, change_range_coeff}, {0.1, 0.05, change_range_coeff * 1.5},
    //        {0.1, 0.05, change_range_coeff}, {0.1, 0.05, change_range_coeff}, {0.1, 0.05, change_range_coeff} };
    //    double intensity = 5;
    //    peakParams = { { 1.9, 110, intensity, 0, 0, 0 },
    //                { 1.44, 162, intensity, 0, 0, 0},
    //                { 2.24, 210, intensity, 0, 0, 0},
    //                { 2.25, 232, intensity, 0, 0, 0},
    //                { 2.10, 242, intensity, 0, 0, 0},
    //                { 2.00, 263, intensity, 0, 0, 0} };
    //}
    else if (type == 200) {
        change_range_coeff = 0.2;
        change_range = { {0.16, 0.1, change_range_coeff}, {0.18, 0.1, change_range_coeff}, {0.4, 0.1, change_range_coeff}, {0.3, 0.1, change_range_coeff}, 
            {0.15, 0.1, change_range_coeff}, {0.2, 0.1, change_range_coeff}, {0.18, 0.1, change_range_coeff}, {0.2, 0.1, change_range_coeff}, {0.18, 0.1, change_range_coeff} };
        index_coff = temp[max_index] - 172;
        if (abs(index_coff) > 50)
            index_coff = copysign(50.0, index_coff);
        peakParams = { { 1.7, temp[max_index] - 79, max_intensity / 10, 0, 0, 0},
                    { 1.35, temp[max_index] - 62, max_intensity / 5, 0, 0, 0},
                    { 1.13, temp[max_index] - 28, max_intensity * 0.72, 0, 0, 0},
                    { 1.01, temp[max_index], max_intensity, 0, 0, 0},
                    { 1.53, temp[max_index] - 5, max_intensity / 6, 0, 0, 0},
                    { 1.00, temp[max_index] + 12, max_intensity * 0.4, 0, 0, 0},
                    { 0.93, temp[max_index] + 55, max_intensity * 0.73, 0, 0, 0},
                    { 1.38, temp[max_index] + 61, max_intensity * 0.27, 0, 0, 0},
                    { 1.41, temp[max_index] + 67, max_intensity * 0.4, 0, 0, 0} };
        //if (out_bounds(peakParams, temp)) {
        //    FOM = 1;
        //    return 1000;
        //}
    }
    else if (type == 300) {
        change_range_coeff = 0.2;
        change_range = { {0.3, 0.03, change_range_coeff}, {0.23, 0.02, change_range_coeff}, {0.4, 0.03, change_range_coeff}, {0.37, 0.03, change_range_coeff}, 
            {0.24, 0.023, change_range_coeff}, {0.2, 0.012, change_range_coeff}, {0.3, 0.024, change_range_coeff}, {0.24, 0.023, change_range_coeff} };
        peakParams = { {1.45, temp[max_index] - 94, max_intensity * 0.072, 0, 0, 0},
                    { 1.17, temp[max_index] - 76, max_intensity * 0.169, 0, 0, 0},
                    { 1.3, temp[max_index] - 55, max_intensity * 0.254, 0, 0, 0},
                    { 1.09, temp[max_index] - 27, max_intensity * 0.176, 0, 0, 0},
                    { 1.31, temp[max_index], max_intensity, 0, 0, 0},
                    { 1.26, temp[max_index] + 28, max_intensity * 0.07, 0, 0, 0},
                    { 1.43, temp[max_index] + 77, max_intensity * 0.26, 0, 0, 0},
                    { 1.31, temp[max_index] + 101, max_intensity * 0.2465, 0, 0, 0} };
        if (out_bounds(peakParams, temp)) {
            FOM = 1;
            return 1000;
        }
    }
    else if (type == 400) {
        change_range_coeff = 0.2;
        change_range = { {0.15, 0.014, change_range_coeff}, {0.21, 0.015, change_range_coeff}, {0.18, 0.013, change_range_coeff} };

        peakParams = { {1.42, temp[max_index] - 27, max_intensity / 1.55, 0, 0, 0},
                    { 1.32, temp[max_index], max_intensity, 0, 0, 0},
                    { 1.36, temp[max_index] + 26, max_intensity / 1.3, 0, 0, 0} };
        //if (out_bounds(peakParams, temp)) {
        //    FOM = 1;
        //    return 1000;
        //}
    }
    else if (type == 700) {
        index_coff = temp[max_index] - 237;
        change_range = { {0.15, 0.03, 0.05}, {0.15, 0.017, 0.05}, {0.14, 0.02, 0.05}, {0.06, 0.013, 0.05} };
        peakParams = { { 1.5, temp[max_index] - 125, max_intensity * 0.0323, 0, 0, 0 },
                    { 1.7, temp[max_index] - 50, max_intensity * 0.0323, 0, 0, 0},
                    { 2.7, temp[max_index], max_intensity, 0, 0, 0},
                    { 5.5, temp[max_index] + 10, max_intensity * 0.0323, 0, 0, 0} };
        if (out_bounds(peakParams, temp)) {
            FOM = 1;
            return 1000;
        }
    }

    else {
        change_range_coeff = 0.2;
        change_range = { {0.12, 0.035, change_range_coeff}, {0.17, 0.016, change_range_coeff}, {0.24, 0.08, change_range_coeff}, {0.14, 0.026, change_range_coeff},
            {0.155, 0.012, change_range_coeff}, {0.1, 0.02, change_range_coeff}, {0.25, 0.053, change_range_coeff} };
        peakParams = { { 1.02, temp[max_index] - 130, max_intensity / 8.55, 0, 0, 0},
                     { 0.96, temp[max_index] - 110, max_intensity / 8.125, 0, 0, 0},
                     { 1.05, temp[max_index] - 93, max_intensity / 5, 0, 0, 0},
                     { 1.40, temp[max_index] - 67, max_intensity / 5.91, 0, 0, 0},
                     { 1.42, temp[max_index] - 53, max_intensity / 13, 0, 0, 0},
                     { 0.96, temp[max_index] - 10, max_intensity / 5, 0, 0, 0},
                     { 0.88, temp[max_index] + 5, max_intensity, 0, 0, 0} };
        if (out_bounds(peakParams, temp)) {
            FOM = 1;
            return 1000;
        }
    }
    orig_peak = peakParams;
    int curveSize = int(curve.size());
    int peakNum = int(peakParams.size());
    
    //find the half max index
    vector<double> temperature = temp;
    //find full width half max points
    find_index(temperature, peakParams);
    //calculate original fom
    double orig_integral = 0.0;
    vector<double> orig_fit(curveSize);
    for (int d = 0; d < curveSize; d++) {
        double fit = 0.0;
        for (int e = 0; e < peakNum; e++) {
            fit += func(temp[d], peakParams[e]);
        }
        orig_integral += fit;
        orig_fit[d] = fit;
    }
    double orig_fom = 0.0;
    for (int f = 0; f < int(curve.size()); ++f) {
        orig_fom += abs(curve[f] - orig_fit[f]) / orig_integral;
    }
    original_fom = orig_fom;
    vector<vector<double>> temp_params = peakParams;
    double k = .000086173303;
    double rate1 = 0.0000001;
    double rate2 = 0.00002;
    double rate3 = 0.00001;
    double current_FOM = orig_fom;
    int iteration = 0;
    //bool matrix to check if should apply change
    vector<vector<bool>> process(peakNum, vector<bool>(3, true));
    //check fom
    bool check = true;
    //check if any still inside allowable range
    bool continue_process = true;
    while (iteration < 500 && check && continue_process && current_FOM > 0) {
        vector<vector<double>> update(peakNum, vector<double>(3, 0.0));
        vector<vector<double>> curves(temp_params.size(), vector<double>(curve.size()));
        vector<double> totalCurve(curve.size());
        for (int i = 0; i < int(curve.size()); ++i) {
            double out;
            double accum = 0.0;
            for (int x = 0; x < int(temp_params.size()); ++x) {
                out = quickFok(temp[i], temp_params[x]);
                curves[x][i] = out;
                accum += out;
            }
            totalCurve[i] = accum;
        }
        //use FWHM to find the left and right half max points
        for (int b = 0; b < peakNum; b++) {
            int TL = temp_params[b][3];
            int TR = temp_params[b][5];
            double energy = temp_params[b][0];
            double Tm = temp_params[b][1];
            double Im = temp_params[b][2];
            for (int index = TL; index < TR + 1; index++) {
                double y = curve[index] - totalCurve[index] + curves[b][index];
                if (y < 0)
                    y = 0;
                double T = temp[index];
                double deriv_E = -2.0 * Im * (-(2.0 * k * T * T * T * exp((energy * (T - Tm)) / (Tm * k * T))) / (energy * energy * Tm * Tm) + (2.0 * Tm * k) /
                    (energy * energy) - (T * (T - Tm) * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm * Tm * k) + (T - Tm) /
                    (Tm * k * T)) * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) /
                    (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) /
                        (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0));
                double deriv_Tm = -2.0 * Im * ((2.0 * T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm * Tm) - (T * T * (1.0 - (2.0 * k * T) / energy)
                    * exp((energy * (T - Tm)) / (Tm * k * T)) * (-(energy * (T - Tm)) / (Tm * Tm * k * T) - energy / (Tm * k * T))) / (Tm * Tm) - (energy * (T - Tm)) / (Tm * Tm * k * T) -
                    energy / (Tm * k * T) - (2.0 * k) / energy) * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) -
                    (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) -
                        (2.0 * Tm * k) / energy + 1.0));
                double deriv_Im = -2.0 * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) / (Tm * Tm) + (energy * (T - Tm)) /
                    (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0) * (y - Im * exp(-(T * T * (1.0 - (2.0 * k * T) / energy) * exp((energy * (T - Tm)) / (Tm * k * T))) /
                    (Tm * Tm) + (energy * (T - Tm)) / (Tm * k * T) - (2.0 * Tm * k) / energy + 1.0));
                update[b][0] += rate1 * deriv_E;
                update[b][1] += rate2 * deriv_Tm;
                update[b][2] += rate3 * deriv_Im;
            }
        }
        //apply changes and check if out range
        for (int c = 0; c < peakNum; c++) {
            if (process[c][0]) {
                temp_params[c][0] -= update[c][0];
                if ((abs(temp_params[c][0] - peakParams[c][0]) / peakParams[c][0]) > change_range[c][0]) {
                    //revert
                    temp_params[c][0] += update[c][0];
                    //stop process it
                    process[c][0] = false;
                }
            }
            if (process[c][1]) {
                temp_params[c][1] -= update[c][1];
                if ((abs(temp_params[c][1] - peakParams[c][1]) / peakParams[c][1]) > change_range[c][1]) {
                    //revert
                    temp_params[c][1] += update[c][1];
                    //stop process it
                    process[c][1] = false;
                }
            }
            if (process[c][2]) {
                temp_params[c][2] -= update[c][2];
                if ((abs(temp_params[c][2] - peakParams[c][2]) / peakParams[c][2]) > change_range[c][2]) {
                    //revert
                    temp_params[c][2] += update[c][2];
                    //stop process it
                    process[c][2] = false;
                }
            }
        }
        //check if all process is false
        continue_process = false;
        for (vector<bool>& b : process) {
            for (int i = 0; i < 3; i++) {
                if (b[i])
                    continue_process = true;
            }
        }
        //calculate new fom
        double new_integral = 0.0;
        vector<double> new_output(curve.size(), 0.0);
        for (int d = 0; d < curveSize; d++) {
            double new_fit = 0.0;
            for (int e = 0; e < peakNum; e++) {
                new_fit += func(temp[d], temp_params[e]);
            }
            new_integral += new_fit;
            new_output[d] = new_fit;
        }
        double new_fom = 0.0;
        for (int f = 0; f < curveSize; ++f) {
            new_fom += abs(curve[f] - new_output[f]) / new_integral;
        }
        if (new_fom > current_FOM)
            check = false;
        else
            current_FOM = new_fom;
        //cout << ".";
        //cout.flush();
        iteration++;
    }
    FOM = current_FOM;
    peakParams = temp_params;
    cout << type << " orig_fom: " << original_fom << " new_fom: " << FOM << endl;
    for (int i = 0; i < static_cast<int>(peakParams.size()); i++) {
        cout << peakParams[i][0] << " " << peakParams[i][1] << " " << peakParams[i][2] << endl;
    }
    return iteration;
    
}

void reformat_name(string filename, string output_dir, pair<vector<double>, vector<double>>& data) {
    //find the heating rate / experiment name
    size_t mid, last;
    last = filename.find_last_of('-');
    if (filename.find_last_of('_') != string::npos)
        mid = filename.find_last_of('_');
    else
        mid = filename.find_last_of('-', last - 1);
    string type = filename.substr(mid + 1, last - mid - 1);
    size_t dot = filename.find_last_of('.');
    string filenumber = filename.substr(last + 1, dot - last);
    int filenum = stoi(filenumber);
    string MLname = "TLD900-" + type + "-" + to_string(filenum) + ".csv";
    ofstream file5;
    string path = output_dir + "/" + MLname;
    file5.open(path);
    file5 << "temp, smoothed";
    file5 << ",\n";
    for (int i = 0; i < int(data.first.size()); i++) {
        file5 << data.first[i] << ",";
        file5 << data.second[i] << ",";
        file5 << "\n";
    }
    file5.close();
}

int main(int argc, char* argv[]) {
    //enable quick mode to run output with machine generated peak detections
    bool output_mode = false;
    if (argc > 1)
        output_mode = (*(argv[1]) == 'q') ? true : false;
    //string to store the input directory path
    string dir;
    //functional string that will be used in the while loop to read in path
    string start = "n";
    //vector to store the processed files, will be used in the statistic function call
    vector<string> filenames;
    //vector of strings that will store all the csv file names in the "dir" directory
    vector<string> files;
    string output_dir = "";
    //sampling rate default to be 0.1
    int time = 0.1;

    //read in the path that contains all the input data
    while (start == "n" || start == "N") {
        cout << "Please enter the full path to directory containing csv formatted emission spectra:" << endl;
        cin >> dir;
        //cast dir in the right format to be opened
        if (dir.back() == '/')
            dir.pop_back();

        //HANDLE FILES
        //call the handle_dir function from FileHandler.cpp, store all csv files in dir to the files vector
        files = handle_dir(dir, output_dir);
        if (!output_mode) {
            cout << "Is this correct and would you like to start processing (y/n)?" << endl;
            cin >> start;
        }
        else
            start = "y";
    }
    //2d vector to store processed data of each input file
    vector<vector<double>> stats(files.size(), vector<double>(0, 0.0));
    int count = 0;

    vector<vector<double>> peakParams;
    vector<vector<vector<double>>> all_peakParam;
    vector<int> samp_rate;
    enum Mode { ALL, EACH, NONE };
    Mode m = Mode::NONE;
    if (!output_mode) {
        //ASK USER PEAK INPUT
        string choice;
        cout << "Do you want to manually input peak locations and heights (y/n)?" << endl;
        cin >> choice;
        if (choice == "y") {
            string repeat;
            cout << "Do you want the same input for all files or different input for each file (all/each)?" << endl;
            cin >> repeat;
            //ask user to input single set of peaks data for all files
            if (repeat == "all") {
                m = Mode::ALL;
                cout << "Please type in data in the format: tmeperature,count,activation energy, press enter for each peak." << endl;
                cout << "Type done when you are finished." << endl;
                peakParams = input_data();
                if (peakParams.empty()) {
                    m = Mode::NONE;
                    cout << "Empty input, switching to automatic peak identification." << endl;
                }
                cout << "Please give the sampling rate for all data.";
                cin >> time;
            }
            else if (repeat == "each") {
                m = Mode::EACH;
                cout << "For each file, please type in peak data in the format: tmeperature,count,activation energy, press enter for each peak." << endl;
                cout << "Type done when you are finished." << endl;
                cout << "Then, please type in the sampling rate." << endl;
                auto i = files.begin();
                for (; i != files.end(); ++i) {
                    string filename = i->substr((i->find_last_of("/\\")) + 1);
                    cout << "Enter peak data for: " << filename << endl;
                    vector<vector<double>> param = input_data();
                    all_peakParam.push_back(param);
                    cout << "Enter sampling rate." << endl;
                    int t;
                    cin >> t;
                    samp_rate.push_back(t);
                }
                if (all_peakParam.empty()) {
                    m = Mode::NONE;
                    cout << "Empty input, switching to automatic peak identification." << endl;
                }
            }
            else {
                cout << "Invalid command!" << endl;
            }
        }
    }
    ofstream file1;
    string path = "C:/Users/jack0/Desktop/report.csv";
    file1.open(path);
    //ofstream output;
    //string path = "C:/Users/jack0/Desktop/graph.txt";
    //output.open(path);
    ofstream file9;
    string outpath = output_dir + "/statistic.csv";
    file9.open(outpath);
    for (int i = 0; i < static_cast<int>(files.size()); i++) {
        //count to see which file is being processed
        int num = 0;
        //make a coopy of peakParam and then customize the copy accoding to the data read in
        vector<vector<double>> peak;
        if (m == Mode::ALL)
            peak = peakParams;
        if (m == Mode::EACH)
            peak = all_peakParam[num];

        //erase the previous temp.csv file to read in new data
        if (files[i].find("temp.csv") != string::npos) {
            files.erase(files.begin() + i);
            continue;
        }
        cout << "----------------------------" << endl << "Processing: ";
        string filename = files[i].substr((files[i].find_last_of("/\\")) + 1);
        cout << filename << " (" << count + 1 << " of " << files.size() << ")" << endl << "Reading in File  .";
        cout.flush();

        //FILE_MANAGER created
        //create a fileManager object that takes in the i/csv path
        File_Manager fileManager = *new File_Manager(files[i], output_dir);
        cout << "." << endl;
        cout.flush();
        //create a pair of two vector data which has first to be temperature data and second to be count data
        pair<vector<double>, vector<double>> data = fileManager.read();
        int length = static_cast<int>(data.second.size());
        int window_size = length * 0.05;
        if (window_size % 2 == 0) {
            window_size += 1;
        }

        //file1 << filename << ",";
        //auto start = chrono::high_resolution_clock::now();

        //Smooth temperature raw data with Savitzky-Golay
        int window = length * 0.05;
        if (window % 2 == 0) {
            window += 1;
        }
        
        vector<double> origin_temp = data.first;
        vector<double> origin_count = data.second;

        vector<double> temp_orig1 = data.first;
        vector<double> temp_orig2 = temp_orig1;

        SG_smooth(temp_orig2, window, 3);

        //for (int i = 0; i < static_cast<int>(temp_orig1.size()); i++) {
        //    data.first[i] = (temp_orig1[i] + temp_orig2[i]) / 2;
        //}
        for (int i = 0; i < static_cast<int>(temp_orig1.size()); i++) {
            data.first[i] = temp_orig2[i];
        }
        
        //REMOVE_SPIKE call
        spike_elim(data.first, data.second, 3, 1.2);
        vector<double> orig_count = data.second;
        //spike_elim(data.first, data.second, 3, 1.05);
        
        //copy two times the count data and run Savitzky-Golay with order 4 and 5, then take the average
        vector<double> orig_count1 = data.second;
        vector<double> orig_count2 = orig_count1;
        //SG_smooth(data.second, window_size, 2);
        SG_smooth(orig_count1, window_size, 4);
        SG_smooth(orig_count2, window_size, 5);
        for (int i = 0; i < static_cast<int>(orig_count1.size()); i++) {
            data.second[i] = (orig_count1[i] + orig_count2[i]) / 2;
        }
        //for (int i = 0; i < int(data.first.size()); i++) {
        //    file1 << origin_temp[i] << ",";
        //    file1 << origin_count[i] << ",";
        //    file1 << data.first[i] << ",";
        //    file1 << orig_count[i] << ",";
        //    file1 << data.second[i] << ",";
        //    file1 << endl;
        //
        //}
        //file1.close();
        

        //remove_tail(data.first, data.second, max_intensity);
        //background_substraction
        //vector<double> t = remove_back(data.first, data.second);
        
        vector<double> smooth_count = data.second;
        double max_intensity = *max_element(data.second.begin(), data.second.end());

        //calculate the curve area by adding the count data
        const double curveArea = accumulate(data.second.begin(), data.second.end(), 0.0);
        //if the curve area is less 2000 then it's not enough for further analysis
        if (curveArea < 2000) {
            files.erase(files.begin() + i);
            i--;
            cout << "file not considered" << endl;
            //remove((dir + "/temp.csv").c_str());
            remove((output_dir + "/temp.csv").c_str());
            continue;
        }
        vector<double> smoothed_count = data.second;


        //output_counts(orig_count, smoothed_count, data.first, filename);
        //testing gradient descent on TLD 100
        //peakparams: activation energy, maxTemp, maxIntensity, TL, TM, TR
        //vector<vector<double>> peakParams;
        //findPeaks(data.first, data.second, peakParams, output_dir);
        //if (int(peakParams.size()) != 4) {
        //    files.erase(files.begin() + i);
        //    i--;
        //    cout << "file not considered" << endl;
        //    //remove((dir + "/temp.csv").c_str());
        //    remove((output_dir + "/temp.csv").c_str());
        //    continue;
        //}
        //cout << "energy, temp, count, TL, TM, TR" << endl;
        //for (int i = 0; i < int(peakParams.size()); i++)
        //    cout << peakParams[i][0] << " " << peakParams[i][1] << " " << peakParams[i][2] << " " << peakParams[i][3] << " " << peakParams[i][4] << " " << peakParams[i][5] << endl;
        //vector<vector<double>> curve;
        //for (int i = 0; i < int(peakParams.size()); ++i) {
        //    curve.push_back(vector<double>(data.first.size(), 0.0));
        //}
        //calculate every temperature's FOK data in each peak fit, accumulate peak areas for each peak in
        //peak_areas and accumulate same temperature's FOK values in all fits to sum

        // uncomment
        vector<vector<double>> peak_param;
        
        //vector<vector<double>> peak_100;
        //vector<vector<double>> orig_peak_100;
        //double fom_100 = -1;
        //double original_fom_100 = 0;
        //int iteration_100 = gd_types(data.first, data.second, peak_100, fom_100, 100, max_intensity, original_fom_100, orig_peak_100);
        //if (isnan(fom_100))
        //    fom_100 = 100;
        //peak_param = peak_100;
        
        
        //vector<vector<double>> peak_101;
        //vector<vector<double>> orig_peak_101;
        //double fom_101 = -1;
        //double original_fom_101 = 0;
        //int iteration_101 = gd_types(data.first, data.second, peak_101, fom_101, 101, max_intensity, original_fom_101, orig_peak_101);

        //vector<vector<double>> peak_200;
        //vector<vector<double>> orig_peak_200;
        //double fom_200 = -1;
        //double original_fom_200 = 0;
        //int iteration_200 = gd_types(data.first, data.second, peak_200, fom_200, 200, max_intensity, original_fom_200, orig_peak_200);
        //
        //
        //if (iteration_200 == 1000) {
        //    files.erase(files.begin() + i);
        //    i--;
        //    cout << "file not considered" << endl;
        //    //remove((dir + "/temp.csv").c_str());
        //    remove((output_dir + "/temp.csv").c_str());
        //    continue;
        //}
        
        //vector<vector<double>> peak_300;
        //vector<vector<double>> orig_peak_300;
        //double fom_300 = -1;
        //double original_fom_300 = 0;
        //int iteration_300 = gd_types(data.first, data.second, peak_300, fom_300, 300, max_intensity, original_fom_300, orig_peak_300);
        
        
        
        //vector<vector<double>> peak_400;
        //double fom_400 = -1;
        //double original_fom_400 = 0;
        //vector<vector<double>> orig_peak_400;
        //int iteration_400 = gd_types(data.first, data.second, peak_400, fom_400, 400, max_intensity, original_fom_400, orig_peak_400);
        
        
        
        //vector<vector<double>> peak_700;
        //vector<vector<double>> orig_peak_700;
        //double fom_700 = -1;
        //double original_fom_700 = 0;
        //int iteration_700 = gd_types(data.first, data.second, peak_700, fom_700, 700, max_intensity, original_fom_700, orig_peak_700);
        
        
        
        vector<vector<double>> peak_900;
        vector<vector<double>> orig_peak_900;
        double fom_900 = -1;
        double original_fom_900 = 0;
        int iteration_900 = gd_types(data.first, data.second, peak_900, fom_900, 900, max_intensity, original_fom_900, orig_peak_900);
        //
        //vector<double> fom_set{ abs(fom_100), abs(fom_200), abs(fom_300), abs(fom_400), abs(fom_700), abs(fom_900) };
        //int min_fom = min_element(fom_set.begin(), fom_set.end()) - fom_set.begin();
        //
        //int adopt = -1;
        //double orig_fom = -1;
        //int iter = -1;
        //double cur_fom = -1;
        //
        //if(min_fom == 0) {
        //    peak_param = peak_100;
        //    adopt = 100;
        //    orig_fom = original_fom_100;
        //    cur_fom = fom_100;
        //    iter = iteration_100;
        //}
        //else if (min_fom == 1) {
        //    peak_param = peak_200;
        //    adopt = 200;
        //    orig_fom = original_fom_200;
        //    cur_fom = fom_200;
        //    iter = iteration_200;
        //}
        //else if (min_fom == 2) {
        //    peak_param = peak_300;
        //    adopt = 300;
        //    orig_fom = original_fom_300;
        //    cur_fom = fom_300;
        //    iter = iteration_300;
        //}
        //else if (min_fom == 3) {
        //    peak_param = peak_400;
        //    adopt = 400;
        //    orig_fom = original_fom_400;
        //    cur_fom = fom_400;
        //    iter = iteration_400;
        //}
        //else if (min_fom == 4) {
        //    peak_param = peak_700;
        //    adopt = 700;
        //    orig_fom = original_fom_700;
        //    cur_fom = fom_700;
        //    iter = iteration_700;
        //}
        //else {
        //    peak_param = peak_900;
        //    adopt = 900;
        //    orig_fom = original_fom_900;
        //    cur_fom = fom_900;
        //    iter = iteration_900;
        //}
        //stop
        peak_param = peak_900;

        //auto stop = chrono::high_resolution_clock::now();
        //auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        //file1 << duration.count() << ",";

        //output_details(output_dir, filename, data.first, data.second, peak_300, orig_peak_300, peak_400, orig_peak_400, peak_900, orig_peak_900, peak_100, orig_peak_100, peak_200, orig_peak_200);

        //auto start1 = chrono::high_resolution_clock::now();
        // LM output
        First_Order_Kinetics FOK_Model = *new First_Order_Kinetics(data, peak_param);
        file9 << filename << ",";
        ////////TLD 700 first temp 0.2
        vector<vector<double>> constrain = {
            {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0, 0},
            {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.05, 0},
            {0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3 }
        }; 
        double res = FOK_Model.glow_curve(file9, constrain);
        if (res >= 1) {
            files.erase(files.begin() + i);
            continue;
        }
        // 
        //output glow peaks
        //vector<vector<double>> returnedPeaks = FOK_Model.return_glow_curve();
        //ofstream file5;
        //string path = output_dir + "/" + filename;
        //file5.open(path);
        //file5 << "temp, smoothed, first, sec, third";
        //file5 << ",\n";
        //for (int i = 0; i < int(data.first.size()); i++) {
        //    file5 << data.first[i] << ",";
        //    file5 << smoothed_count[i] << ",";
        //    for (int j = 0; j < int(peak_param.size()); j++) {
        //        file5 << returnedPeaks[j][i] << ",";
        //    }
        //    file5 << "\n";
        //
        //}
        //file5.close();
        //auto stop1 = chrono::high_resolution_clock::now();
        //auto duration1 = chrono::duration_cast<chrono::milliseconds>(stop1 - start1);
        //if (res == -1) {
        //    file1 << duration1.count() << "," << "failed" << endl;
        //}
        //else {
        //    file1 << duration1.count() << "," << endl;
        //}
        //
        //for(int i = 0; i < static_cast<int>(peak_param.size()); i++) {
        //    for (int j = 0; j < static_cast<int>(peak_param[0].size()); j++) {
        //        cout << peak_param[i][j] << " ";
        //    }
        //    cout << endl;
        //}
        count++;
    }
    //file1.close();
    return 0;
}

