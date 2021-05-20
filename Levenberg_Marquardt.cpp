//
//  First_Order_kinetics.cpp
//  GlowCurveAnalsys
//
//  Initially Created by jeremy hepker on 1/27/19.
//
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//

#include "Levenberg_Marquardt.hpp"

using namespace std;

First_Order_Kinetics::First_Order_Kinetics(std::pair<std::vector<double>,std::vector<double>> data, std::vector<std::vector<double>> peakParams)
    :count_data(data.second),temp_data(data.first), peakParams(peakParams){
    //calculate 1st derivative and populate orig_sig_deriv
    deriv(temp_data, count_data, orig_sig_deriv);
    };

/*---------------------------------Main Deconvolution---------------------------------------*/
//calculate FOM and the area under each curve fit
double First_Order_Kinetics::glow_curve(ofstream& file, vector<vector<double>>& constrain){
    double FOM = 1.0;
    vector<double> sum(count_data.size(), 0.0);
    double integral = 0.0;
    cout << ".";
    cout.flush();
    //call LevenbergMarquardt from this file
    LevenbergMarquardt(count_data, peakParams, FOM, constrain);
    //gradient_Descent(count_data, peakParams, FOM);
    //LevenbergMarquardt(count_data, peakParams, FOM);
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
    //for(int i = 0; i < int(peak_areas.size());++i ){
    //    if(peak_areas[i] < 200.0){
    //        peak_areas.erase(peak_areas.begin()+i);
    //        peakParams.erase(peakParams.begin()+i);
    //        --i;
    //    }
    //}
    //glow_curves.push_back(sum);
    //output the area under each peak
    for(int i = 0; i < int(peakParams.size()); ++i){
        cout << "----- Area Under Curve #" << i+1 << " :" << peak_areas[i] << " -----" << endl;
    }
    curve_areas = peak_areas;
    totalArea = integral;
    //file << FOM << ",";
    //for (vector<double> v : peakParams) {
    //    file << v[0] << "," << v[1] << "," << v[2] << ",";
    //}
    //for (double d : peak_areas) {
    //    file << d << ",";
    //}
    //file << endl;
    return FOM;
};

/*-------------------------First Order Kinetics Function--------------------------------*/

//this is the same as FOKModel.cpp, calculate FOK data
double First_Order_Kinetics::Func2(const double input, const vector<double> params){
    double T=0.0;
    double I_t = 0.0;
    double energy = params[0];
    double Tm = params[1]+273.15;
    double dm = (2.0*k*(Tm))/ energy;
    double Im = params[2];
    T = double(input+273.15);
    I_t = Im*exp(1.0 +(energy/(k*T))*((T-Tm)/Tm)-((T*T)/(Tm*Tm))*exp((energy/(k*T))*((T-Tm)/Tm))*(1.0-((2.0*k*T)/energy))-dm);
    return I_t;
}

//-------------------------Levenburg Marquardt METHOD----------------------------------------//
// use Levenberg-Marquardt method to further fit the curve and calculate figure of merit for the curve fitting
void First_Order_Kinetics::LevenbergMarquardt(const vector<double> &curve, vector<vector<double>> &params, double &FOM, vector<vector<double>>& constrain){
    auto start = chrono::high_resolution_clock::now();
    //create singlePeak vector that is same size as curve
    vector<double> singlePeak(curve.size(), 0.0);
    int curveSize = int(curve.size());
    int peakSize = int(params.size());
    int d = 1;
    int main_hold = 0;
    double main_FOM = FOM;
    vector<double> orig_energy(params.size());
    vector<double> orig_temp(params.size());
    vector<double> orig_height(params.size());
    for (int i = 0; i < int(params.size()); i++) {
        orig_energy[i] = params[i][0];
    }
    for (int i = 0; i < int(params.size()); i++) {
        orig_temp[i] = params[i][1];
    }
    for (int i = 0; i < int(params.size()); i++) {
        orig_height[i] = params[i][2];
    }
    while(FOM > .01){
        main_FOM = FOM;
        if(main_hold > 3){
            break;
        }
        //param_num is used to indicate which value (energy, temp, count) will be modifed in the L-M process
        for(int param_num = 0; param_num < 3 ; ++param_num){
            vector<double> temp_params;
            vector<double> temp_output(curve.size(), 0.0);
            //push ith value in params to temp_params
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
            //2d vector to store each point's derivative under each peak curve
            vector<vector<double>> Jf_T(peakSize, vector<double>(curveSize,0.0));
            vector<double> error(curveSize,0.0);
            vector<vector<double>> H;
            double lambda = 0.01;
            double updateJ = 1;
            //double e = 0.0;
            int i = 0;
            int inner_hold = 0;
            //the FOM gets better or it stops at 300 iteration
            while(FOM > .02 && i < 300){
                auto stop = chrono::high_resolution_clock::now();
                auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
                if (duration.count() > 240000) {
                    main_hold = 4;
                    break;
                }
                if(updateJ == 1){
                    //Evaluate the jacobian matrix at the current paramater.
                    vector<double>process_parms(3, 0.0);
                    double integral = 0.0;
                    for(int j = 0; j < curveSize; j++) {
                        double output = 0.0;
                        for(int k = 0; k < peakSize; ++k){
                            process_parms[param_num] = temp_params[k];
                            process_parms[other_param1] = params[k][other_param1];
                            process_parms[other_param2] = params[k][other_param2];
                            Jf_T[k][j] = Deriv2(temp_data[j], process_parms, param_num);
                            //output accumulates the single point's fitted counts
                            output += Func2(temp_data[j], process_parms);
                        }
                        temp_output[j] = output;
                        //error stores the difference in orginal count and the accumulated fitted count
                        error[j] = curve[j] - output;
                        integral += output;
                    }
                    //this is where the figure of merit is calculated before each step
                    FOM = 0.0;
                    for(int z = 0; z < int(curve.size()); ++z){
                        FOM += abs(curve[z] - temp_output[z])/integral;
                    }
                    //Calculate the hessian matrix
                    vector<vector<double>> Jf(Jf_T[0].size(), vector<double>(Jf_T.size(),0.0));
                    //Jf is transpose of Jf_T
                    transpose(Jf_T, Jf, int(Jf_T.size()), int(Jf_T[0].size()));
                    //multiplication of the jacobian matrix and its transpose
                    H = multiply(Jf_T,Jf);
                    //e = dotProduct(error, error);
                }

                //apply the damping factor to the hessian matrix
                //I is a diagonal matrix with diagonal equals lambda which is the damping vector
                vector<vector<double>> I = Identity(peakSize, lambda);
                vector<vector<double>> H_lm(peakSize, vector<double>(peakSize,0.0));
                for(int j = 0; j < int(H.size()); ++j){
                    for(int s = 0; s < int(H.size()); ++s){
                        H_lm[j][s] = H[j][s] + I[j][s];
                    }
                }
                invert(H_lm, true);
                //multiply jacobian matrix with error matrix
                vector<double> Jf_error = vec_matrix_multi(Jf_T, error);
                vector<double> delta = vec_matrix_multi(H_lm, Jf_error);
                vector<double> t_params = temp_params;
                //update the change to the original data
                for(int x = 0; x < int(delta.size()); ++x){
                    t_params[x] += delta[x];
                    //if (param_num == 0 && (abs((t_params[x] - orig_energy[x]) / orig_energy[x]) > 0.17)) {
                    //    t_params[x] -= delta[x];
                    //}
                    //if (param_num == 0 && x == 1 && (abs((t_params[x] - orig_energy[x]) / orig_energy[x]) > 0.5)) {
                    //    t_params[x] -= delta[x];
                    //}
                    //if (param_num == 2 && (abs((t_params[x] - orig_height[x]) / orig_height[x]) > 0.4)) {
                    //    t_params[x] -= delta[x];
                    //}

                    //if(param_num == 0 && constrain[0][x] != 0 && (abs((t_params[x] - orig_energy[x]) / orig_energy[x]) > constrain[0][x]))
                    //    t_params[x] -= delta[x];
                    //if (param_num == 1 && constrain[1][x] != 0 && (abs((t_params[x] - orig_temp[x]) / orig_temp[x]) > constrain[1][x]))
                    //    t_params[x] -= delta[x];
                    //if (param_num == 2 && constrain[2][x] != 0 && (abs((t_params[x] - orig_height[x]) / orig_height[x]) > constrain[2][x]))
                    //    t_params[x] -= delta[x];


                    //if (param_num == 0 && x == 1) {
                    //    if (abs((t_params[x] - orig_energy[x]) / orig_energy[x]) > 0.02) {
                    //        t_params[x] -= delta[x];
                    //    }
                    //}
                    //if (param_num == 0 && x == 2) {
                    //    if (abs((t_params[x] - orig_energy[x]) / orig_energy[x]) > 0.06) {
                    //        t_params[x] -= delta[x];
                    //    }
                    //}
                    //if (param_num == 2) {
                    //    if (x == 3) {
                    //        if (abs((t_params[x] - orig_height[x]) / orig_height[x]) > 0.15) {
                    //            t_params[x] -= delta[x];
                    //        }
                    //    }
                    //    else if (x == 2 && ((t_params[x] - orig_height[x]) / orig_height[x]) < -0.1) {
                    //        t_params[x] -= delta[x];
                    //    }
                    //}
                }
                double integral = 0.0;
                //Evaluate the total distance error at the updated paramaters.
                vector<double> temp_error(curveSize,0.0);
                vector<double> t_param(3,0.0);
                //for every point under every peak, with the new parameter recalculate count data
                for(int j = 0; j < curveSize; j++){
                    double output = 0.0;
                    for(int k = 0; k < peakSize;++k){
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
                //calculate figure of merit after the step
                double temp_FOM = 0.0;
                if (integral == 0) {
                    temp_FOM = 1;
                    FOM = 1;
                    break;
                }
                else {
                    for (int z = 0; z < int(curve.size()); ++z) {
                        temp_FOM += abs(curve[z] - temp_output[z]) / integral;
                    }
                }
                //if new fom is smaller then proceed to calculate a new jacobian with a smaller damping parameter
                if(temp_FOM < FOM){
                    FOM = temp_FOM;
                    lambda /= 10;
                    temp_params = t_params;
                    updateJ = 1;
                    inner_hold = 0;
                //if new fom is not better, then continue with the current jacobian matrix and use a larger damping parameter
                }else{
                    inner_hold += 1;
                    updateJ = 0;
                    lambda *= 10;
                }
                if(inner_hold > 25) 
                    i = 500;
                ++i;
            }
            //update params to contain new fitted data in temp_params
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

//populate decon_sig_deriv for the sum of the deconvolute curve
void First_Order_Kinetics::update_deriv() {
    vector<double> dir = glow_curves.back();
    deriv(temp_data, dir, decon_sig_deriv);
}

//calculate figure of merit from derivative of original count and derivative of accumulated count after curve fitting
void First_Order_Kinetics::deriv_FOM() {
    double total = 0.0;
    double fom = 0.0;
    //calculate total derivative of original count derative
    for (double i : orig_sig_deriv) {
        total += i;
    }
    for (int i = 0; i < static_cast<int>(orig_sig_deriv.size()); i++) {
        fom += abs(orig_sig_deriv[i] - decon_sig_deriv[i]) / total;
    }
    fom *= 100;
}

void First_Order_Kinetics::gradient_Descent(const vector<double>& curve, vector<vector<double>>& peakParams, double& FOM) {
    //temperary vector to store peak data
    vector<vector<double>> temp_params = peakParams;
    //temperary vector to store accumulated fitted count
    vector<double> temp_output(curve.size(), 0.0);
    int curveSize = int(curve.size());
    int peakNum = int(peakParams.size());
    double rate1 = 0.00000001;
    double rate2 = 0.00002;
    double rate3 = 0.001;
    double current_FOM = FOM;
    int iteration = 0;
    int main_hold = 0;
    while (iteration < 300 && main_hold < 3) {
        //calculate initial FOM
        //double integral = 0.0;
        //for each point calculated the accumulated fitted count using FOK model
        //for (int i = 0; i < curveSize; i++) {
        //    double fit = 0.0;
        //    for (int j = 0; j < peakNum; j++) {
        //        fit += Func2(temp_data[i], temp_params[j]);
        //    }
        //    integral += fit;
        //    temp_output[i] = fit;
        //}
        //FOM = 0.0;
        //for (int k = 0; k < curveSize; ++k) {
        //    FOM += abs(curve[k] - temp_output[k]) / integral;
        //}
        //use gradient descent to calculate change in peak parameter
        //vector<vector<double>> temp_params2 = temp_params;
        vector<vector<double>> update(peakNum, vector<double>(3, 0.0));
        //use FWHM to find the left and right half max points
        find_index(temp_data, temp_params);
        for (int b = 0; b < peakNum; b++) {
            int TL = temp_params[b][3];
            int TR = temp_params[b][5];
            double energy = temp_params[b][0];
            double Tm = temp_params[b][1];
            double Im = temp_params[b][2];
            for (int index = TL; index < TR + 1; index++) {
                double y = curve[index];
                double T = temp_data[index];
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
        for(int i = 0; i < peakNum; i++) {
            if (update[i][0] < 0.001 || update[i][1] < 0.001 || update[i][2] < 0.001)
                break;
        }
        //apply the change to the peak paramter
        for (int c = 0; c < peakNum; c++) {
            temp_params[c][0] -= update[c][0];
            temp_params[c][1] -= update[c][1];
            temp_params[c][2] -= update[c][2];
        }
        for (int j = 0; j < peakNum; j++) {
            cout << update[j][0] << ", " << update[j][1] << ", " << update[j][2] << endl;
            cout << temp_params[j][0] << ", " << temp_params[j][1] << ", " << temp_params[j][2] << endl;
        }
        //re-calculate the new FOM
        //double new_integral = 0.0;
        //vector<double> new_output(curve.size(), 0.0);
        //for (int d = 0; d < curveSize; d++) {
        //    double new_fit = 0.0;
        //    for (int e = 0; e < peakNum; e++) {
        //        new_fit += Func2(temp_data[d], temp_params2[e]);
        //    }
        //    new_integral += new_fit;
        //    new_output[d] = new_fit;
        //}
        //current_FOM = 0.0;
        //for (int f = 0; f < curveSize; ++f) {
        //    current_FOM += abs(curve[f] - new_output[f]) / new_integral;
        //}
        //if (current_FOM < FOM) {
        //    //if no significant improvement then only iterate 3 times
        //    if (abs(current_FOM - FOM) < (1e-4)) {
        //        main_hold += 1;
        //    }
        //    FOM = current_FOM;
        //    temp_params = temp_params2;
        //    iteration++;
        //    cout << ".";
        //    cout.flush();
        //}
        //else {
        //    iteration++;
        //    cout << ".";
        //    cout.flush();
        //}
        cout << ".";
        cout.flush();
        iteration++;
    }
    double new_integral = 0.0;
    vector<double> new_output(curve.size(), 0.0);
    for (int d = 0; d < curveSize; d++) {
        double new_fit = 0.0;
        for (int e = 0; e < peakNum; e++) {
            new_fit += Func2(temp_data[d], temp_params[e]);
        }
        new_integral += new_fit;
        new_output[d] = new_fit;
    }
    current_FOM = 0.0;
    for (int f = 0; f < curveSize; ++f) {
        current_FOM += abs(curve[f] - new_output[f]) / new_integral;
    }
    FOM = current_FOM;
    swap(peakParams, temp_params);
}