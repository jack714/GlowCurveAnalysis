//
//  First_Order_kinetics.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 1/27/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//

#include "Levenberg–Marquardt.hpp"
//#include "smartPeakDetect.hpp"

using namespace std;

First_Order_Kinetics::First_Order_Kinetics(std::pair<std::vector<double>,std::vector<double>> data, std::vector<std::vector<double>> peakParams)
    :count_data(data.second),temp_data(data.first), peakParams(peakParams){
       /*
       populate vector called orig_sig_deriv, derivative from data using central difference/5 point stencil
       */
    };

/*---------------------------------Main Deconvolution---------------------------------------*/
//calculate FOM and the area under each curve fit
double First_Order_Kinetics::glow_curve(){
    double FOM = 1.0;
    vector<double> sum(count_data.size(), 0.0);
    double integral = 0.0;
    cout << ".";
    cout.flush();
    //call LevenbergMarquardt from this file
    LevenbergMarquardt(count_data, peakParams, FOM);
    cout<<".";
    cout.flush();
    if(FOM > 1.0){
        cout << "." << endl << "----- Levenberg-Marquardt failed to converged -----" << endl;
        cout << "data most likely contains to much noise, or is improperly formatted" << endl;
        return -1;
    }
    cout << "." << endl << "----- Levenberg-Marquardt converged to a FOM of " << (FOM*100) << "% -----" <<endl;
    vector<double> peak_areas = vector<double>(peakParams.size(), 0.0);
    //set glow_curves vector to have the same number of vectors as peakParams's size which is number of peaks
    for(int i = 0;i < int(peakParams.size()); ++i){
        glow_curves.push_back(vector<double>(temp_data.size(),0.0));
    }
    //calculate every temperature's FOK data in each peak fit, accumulate peak areas for each peak in
    //peak_areas and accumulate same temperature's FOK values in all fits to sum
    for(int i = 0; i < int(temp_data.size()); ++i ){
        double output = 0.0;
        for(int x = 0; x < int(peakParams.size()); ++x){
            double out = Func2(temp_data[i], peakParams[x]);
            peak_areas[x] += out;
            output += out;
            //glow_curves contains FOK value for each temperature under each peak's fitting
            glow_curves[x][i] = out;
        }
        sum[i] = output;
        integral += output;
    }
    //if the area under a peak is less than 200 then discard that peak data
    for(int i = 0; i < int(peak_areas.size());++i ){
        if(peak_areas[i] < 200.0){
            peak_areas.erase(peak_areas.begin()+i);
            peakParams.erase(peakParams.begin()+i);
            --i;
        }
    }
    glow_curves.push_back(sum);
    //output the area under each peak
    for(int i = 0; i < int(peakParams.size()); ++i){
        cout << "----- Area Under Curve #" << i+1 << " :" << peak_areas[i] << " -----" << endl;
    }
    curve_areas = peak_areas;
    totalArea = integral;
    return FOM;
};

/*-------------------------First Order Kinetics Function--------------------------------*/

//this is the same as FOKModel.cpp, calculate FOK data
double First_Order_Kinetics::Func2(const double input, const vector<double> params){
    double T=0.0;
    double I_t = 0.0;
    double E = params[0];
    double Tm = params[1]+273.15;
    double dm = (2.0*k*(Tm))/E;
    double Im = params[2];
    T = double(input+273.15);
    I_t = Im*exp(1.0 +(E/(k*T))*((T-Tm)/Tm)-((T*T)/(Tm*Tm))*exp((E/(k*T))*((T-Tm)/Tm))*(1.0-((2.0*k*T)/E))-dm);
    return I_t;
}

//-------------------------Levenburg Marquardt METHOD----------------------------------------//
// use Levenberg-Marquardt method to further fit the curve and calculate figure of merit for the curve fitting
void First_Order_Kinetics::LevenbergMarquardt(const vector<double> &curve, vector<vector<double>> &params, double &FOM){
    //create singlePeak vector that is same size as curve
    vector<double> singlePeak(curve.size(), 0.0);
    int curveSize = int(curve.size());
    int curentCurve = int(params.size());
    int d = 1;
    int main_hold = 0;
    double main_FOM = FOM;
    while(FOM > .01){
        main_FOM = FOM;
        if(main_hold > 3){
            break;
        }
        for(int param_num = 0; param_num < 3 ; ++param_num){
            vector<double> temp_params;
            vector<double> temp_output(curve.size(), 0.0);
            //push everything in params to temp_params
            for(int i = 0; i < int(params.size());++i){
                temp_params.push_back(params[i][param_num]);
            }
            int other_param1 = 0;
            int other_param2 = 0;
            if(param_num == 0){
                other_param1 = 1;
                other_param2 = 2;
            }else if(param_num == 1){
                other_param1 = 0;
                other_param2 = 2;
            }else{
                other_param1 = 0;
                other_param2 = 1;
            }
            vector<vector<double>> Jf_T(curentCurve, vector<double>(curveSize,0.0));
            vector<double> error(curveSize,0.0);
            vector<vector<double>> H;
            double lambda = 0.01;
            double updateJ = 1;
            //double e = 0.0;
            int i = 0;
            int inner_hold = 0;
            while(FOM > .02 && i < 300){
                if(updateJ == 1){
                    //Evaluate the jacobian matrix at the current paramater.
                    vector<double>t_parms(3, 0.0);
                    double integral = 0.0;
                    for(int j=0; j < curveSize; j++) {
                        double output = 0.0;
                        for(int k = 0; k < curentCurve; ++k){
                            t_parms[param_num] = temp_params[k];
                            t_parms[other_param1] = params[k][other_param1];
                            t_parms[other_param2] = params[k][other_param2];
                            Jf_T[k][j] = Deriv2(temp_data[j],t_parms, param_num);
                            output += Func2(temp_data[j],t_parms);
                        }
                        temp_output[j] = output;
                        error[j] = curve[j] - output;
                        integral += output;
                    }
                    //this is where the figure of merit is calculated before each step
                    FOM = 0.0;
                    for(int z = 0; z < int(curve.size()); ++z){
                        FOM += abs(curve[z] - temp_output[z])/integral;
                    }
                    //Calculate the hessian mat
                    vector<vector<double>> Jf(Jf_T[0].size(), vector<double>(Jf_T.size(),0.0));
                    transpose(Jf_T, Jf, int(Jf_T.size()), int(Jf_T[0].size()));
                    H = multiply(Jf_T,Jf);
                    //e = dotProduct(error, error);
                }

                //apply the damping factor to the hessian matrix
                vector<vector<double>> I = Identity(curentCurve, lambda);
                vector<vector<double>> H_lm(curentCurve, vector<double>(curentCurve,0.0));
                for(int j = 0; j < int(H.size()); ++j){
                    for(int s = 0; s < int(H.size()); ++s){
                        H_lm[j][s] = H[j][s] + I[j][s];
                    }
                }
                invert(H_lm, true);
                vector<double> Jf_error = vec_matrix_multi(Jf_T,error);
                vector<double> delta = vec_matrix_multi(H_lm,Jf_error);
                vector<double> t_params = temp_params;
                for(int x = 0; x < int(delta.size()); ++x){
                    t_params[x] += delta[x];
                }
                double integral = 0.0;
                //Evaluate the total distance error at the updated paramaters.
                vector<double> temp_error(curveSize,0.0);
                vector<double> t_param(3,0.0);
                for(int j=0; j < curveSize; j++){
                    double output = 0.0;
                    for(int k = 0; k < curentCurve;++k){
                        t_param[param_num] = t_params[k];
                        t_param[other_param1] = params[k][other_param1];
                        t_param[other_param2] = params[k][other_param2];
                        double peak = Func2(temp_data[j],t_param);
                        output += peak;
                    }
                    temp_output[j] = output;
                    integral += output;
                    temp_error[j] = curve[j] - output;
                }
                //glow_curves.push_back(temp_output);
                //calculate figure of merit after the step
                double temp_FOM = 0.0;
                for(int z = 0; z < int(curve.size()); ++z){
                    temp_FOM += abs(curve[z] - temp_output[z])/integral;
                }
                //double temp_e = dotProduct(temp_error, temp_error);
                if(temp_FOM < FOM){
                    lambda /= 10;
                    temp_params = t_params;
                    //e = temp_e;
                    updateJ = 1;
                    inner_hold = 0;
                }else{
                    inner_hold += 1;
                    updateJ = 0;
                    lambda *= 10;
                }
                if(inner_hold > 25) i = 500;
                ++i;
            }
            for(int i = 0; i < int(temp_params.size());++i){
                params[i][param_num]= temp_params[i];
            }
        }
        if(abs(main_FOM - FOM) < (1e-4)){
            main_hold += 1;
        }else{
            main_hold = 0;
        }
        ++d;
        cout<<".";
        cout.flush();
    }
}
