//
//  smartPeakDetect.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 7/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "smartPeakDetect.hpp"

// Sudo Main function for smartPeakDetect, take in temperature and count data then record peaks in peakParams
void findPeaks( std::vector<double>& x, std::vector<double>& y,
                std::vector<std::vector<double>>& peakParams, std::string output_dir )
{
    std::vector<double> xNew = x, yNew= y;
    std::vector<std::vector<double>> peaks;
    std::vector<int> maximums, minimum, inflections;
    std::vector<double> firstDir( x.size( ), 0.0 );
    std::vector<double> secDir( x.size( ), 0.0 );
    //call firstDeriv function from this file and populate firstDir vector with first derivative data
    firstDeriv( xNew, yNew, firstDir );
    //call secDeriv function from this file and populate secDir vector with second derivative data
    secDeriv( xNew, yNew, secDir );
    //call smartPoints from this file and populate maximum, minimum, inflections, further process the peaks in maximum
    smartPoints( xNew, yNew, minimum, maximums,firstDir, secDir, inflections );
    //call pointsParams from this file and populate peakParams with temperature, half left, middle, right max index
    pointsParams( xNew, yNew, maximums, minimum, peakParams );
    // assign peakParams[i][2] to be half max middle count
    for( int i = 0 ; i < int( peakParams.size( ) ) ; i++ )
    {
        peakParams[i][2] = yNew[peakParams[i][4]];
    }
    nonMaxPeaks( xNew, yNew, secDir, maximums, minimum, peakParams, output_dir );
    
    minimum.clear( );
    for( int j = 0 ; j < int( peakParams.size( ) ) ; j++ )
    {
        minimum.push_back( peakParams[j][3] );
        minimum.push_back( peakParams[j][5] );
    }
    printFindings( xNew, yNew, minimum, maximums, inflections, output_dir);
    for( auto i = peakParams.begin( ) ; i != peakParams.end( ) ; i++ )
    {
        std::vector<double> peak( x.size( ), 0.0 );
        FOKModel( xNew, peak, i->at( 1 ), i->at( 2 ), i->at( 0 ) );
        peaks.push_back( peak );
    }
    write( peaks, yNew, xNew, output_dir);
}

// Calculate First Derivatives with temperature and data input
void firstDeriv( std::vector<double>& x, std::vector<double>& y,
                 std::vector<double>& derivative )
{
    int size = int( x.size( ) ) - 2;
    // Five Point Stencil Method
    for( int i = 2 ; i < size ; i++ )
    {
        derivative[i-2] = ( -y[i+2] + 8.0 * y[i+1] - 8.0 * y[i-1] + y[i-2] ) / 12.0;
    }
    size = int( derivative.size( ) ) - 1;
    int sign = 4, lastSign = 0;
    bool positive = true;
    for( int i = 0 ; i < size ; i++ )
    {
        if( derivative[i] < 0.0 && positive && sign >= 4 )
        {
            positive = false;
            lastSign = sign;
            sign = 1;
        }
        else if( derivative[i] < 0.0 && positive && sign < 4 )
        {
            for( int j = 1 ; j <= sign ; j++ )
            {
                derivative[i - j] = -( derivative[i - j] );
                positive = false;
            }
            sign = lastSign + sign;
        }
        else if( derivative[i] > 0.0 && !positive && sign < 4 )
        {
            for( int j = 1 ; j <= sign ; j++ )
            {
                derivative[i - j] = abs( derivative[i - j] );
                positive = true;
            }
            sign = lastSign + sign;
        }
        else if( derivative[i] > 0.0 && !positive && sign >= 4 )
        {
            positive = true;
            lastSign = sign;
            sign = 1;
        }
        else
        {
            sign++;
        }
    }
}

// Calculate Second Derivatives with temperature and data input
void secDeriv( std::vector<double>& x, std::vector<double>& y,
               std::vector<double>& derivative )
{
    int size = int( x.size( ) ) - 2;
    // Five Point Stencil Method
    for( int i = 2 ; i < size; i++ )
    {
        derivative[i-2] = ( -y[i+2] + 16.0 * y[i+1] - 30.0 *
                            y[i] + 16.0 * y[i-1] - y[i-2] ) / 12.0;
    }
    size = int( derivative.size( ) ) - 1;
    int sign = 3, lastSign = 0;
    bool positive = true;
    for( int i = 0 ; i < size ; i++ )
    {
        if( derivative[i] < 0.0 && positive && sign >= 3 )
        {
            positive = false;
            lastSign = sign;
            sign = 1;
        }
        else if( derivative[i] < 0.0 && positive && sign < 3 )
        {
            for( int j = 1 ; j <= sign ; j++ )
            {
                derivative[i - j] = -( derivative[i - j] );
                positive = false;
            }
            sign = lastSign + sign;
        }
        else if( derivative[i] > 0.0 && !positive && sign < 3 )
        {
            for( int j = 1 ; j <= sign ; j++ )
            {
                derivative[i - j] = abs( derivative[i - j] );
                positive = true;
            }
            sign = lastSign + sign;
        }
        else if( derivative[i] > 0.0 && !positive && sign >= 3 )
        {
            positive = true;
            lastSign = sign;
            sign = 1;
        }
        else
        {
            sign++;
        }
    }
}

// Curve Fitting and Vetting, populate minimum, maxima, and inflecPnt vectors
void smartPoints( std::vector<double>& x,
                  std::vector<double>& y,
                  std::vector<int>& minimum,
                  std::vector<int>& maxima,
                  std::vector<double> derivative,
                  std::vector<double> secDerivative,
                  std::vector<int>& inflectPnt )
{
    int size = int( derivative.size( ) );
    minimum.push_back( 0 );
    
    // Pushes all current Maxes and Mins indexes into vector
    for( int i = 1 ; i < size ; i++ )
    {
        //maxima is recorded as first derivative change from positive to negative
        if( derivative[i] < 0.0 && derivative[i - 1] > 0.0 )
        {
            maxima.push_back( i );
        }
        //minimum is recorded as first derivative change from negative to positive
        if( derivative[i] > 0.0 && derivative[i - 1] < 0.0 )
        {
            minimum.push_back( i );
        }
    }
    minimum.push_back( int( x.size( ) ) - 1 );
    
    size = int( secDerivative.size( ) );
    // Pushes inflection points into vector, where second derivative changes sign
    for( int i = 1 ; i < size ; i++ )
    {
        if( secDerivative[i] < 0.0 && secDerivative[i - 1] > 0.0 )
        {
            inflectPnt.push_back( i );
        }
        if( secDerivative[i] > 0.0 && secDerivative[i - 1] < 0.0 )
        {
            inflectPnt.push_back( i );
        }
    }
    
    // Original Value = 10
    int adj_peak_distance = 10;
    // Original Value = 5
    int near_point_height = 5;
    // Original Value = 10
    //int low_threshold = 10;
    
    // loop through maximas and refine peaks
    for( int i = 1 ; i < int( maxima.size( ) ) - 1 ; i++ )
    {
        // if peak is not as tall as adjacent peaks then remove that peak
        if( y[maxima[i]] < y[maxima[i + 1]] && y[maxima[i]] < y[maxima[i - 1]] )
        {
            maxima.erase( maxima.begin( ) + i );
            i--;
            continue;
        }
        
        // if left adjacent peak is not far enough from current peak
        if( abs( x[maxima[i]] - x[maxima[i - 1]] ) <= adj_peak_distance )
        {
            //if right adjacent peak is not far enough from current peak then remove peaks before and after
            if( abs( x[maxima[i]] - x[maxima[i + 1]] ) <= adj_peak_distance )
            {
                maxima.erase( maxima.begin( ) + i + 1 );
                maxima.erase( maxima.begin( ) + i - 1 );
                i--;
                continue;
            }
            //if right adjacent peak is far enough then remove the current peak
            else
            {
                maxima.erase( maxima.begin( ) + i );
                i--;
                continue;
            } // else
        } // if
    } // for
    
    // loop through all peaks and further check height
    for( int i = 0 ; i < int( maxima.size( ) ) ; i++ )
    {
        // if a peak is not as tall as a previous peak or is less than a certain threshold
        //|| y[maxima[i]] < low_threshold
        if( ( y[maxima[i]] < y[maxima[i] - near_point_height] ))
        {
            maxima.erase( maxima.begin( ) + i );
            i--;
            continue;
        }
    }
}

// Poplate peakParams with activation data, temperature, half max left, middle, and right index
void pointsParams( std::vector<double>& x,
                   std::vector<double>& y,
                   std::vector<int>& maxima,
                   std::vector<int>& minima,
                   std::vector<std::vector<double>>& peakParams )
{
    int TM_index = 0, TR_index = 0, TL_index = 0;
    double half_intensity = 0;
    
    // Find the range of the curve
    for( auto i = maxima.begin( ) ; i != maxima.end( ) ; i++ )
    {
        int minLeft = 0, minRight = int( x.size( ) );
        // Look for left adjacent min
        for( int k = *i ; k > 0 ; k-- )
        {
            //iterate through minima and try to find min to the left
            for( auto m = minima.begin( ) ; m != minima.end( ) ; m++ )
            {
                if( k == *m )
                {
                    minLeft = k;
                    break;
                }
            }
            // breaks if a min was found
            if ( minLeft != 0 )
            {
                break;
            }
        }
        // Look for right adjacent min
        for( int k = *i ; k < int( x.size( ) ) ; k++ )
        {
            //iterate through minima and try to find min to the right
            for( auto m = minima.begin( ) ; m != minima.end( ) ; m++ )
            {
                if( k == *m )
                {
                    minRight = k;
                    break;
                }
            }
            if ( minRight != int( x.size( ) ) )
            {
                break;
            }
        }
        peakParams.push_back( std::vector<double>( 6, 0.0 ) );
        auto TL = y.begin( );
        auto TR = y.begin( );
        auto TM = TL + *i;
        half_intensity = *TM / 2.0;
        // Finds half max points
        TL = lower_bound( TL + minLeft, TM, half_intensity, std::less<double>( ) );
        TR = lower_bound( TM, TR + minRight, half_intensity, std::greater<double>( ) );
        // Calculates half width half max
        int diff1 = int( TM - TL );
        int diff2 = int( TR - TM );
        if( TR == y.end( ) )
        {
            diff2 += diff1;
        }
        TM_index = int( TM - y.begin( ) );
        
        // Centering, finds smaller out of left and right side and adjust accordingly
        if( diff1 <= diff2 )
        {
            TL_index = int( TL - y.begin( ) );
            TR_index = int( TM_index + ( TM_index - TL_index ) );
        }
        else
        {
            TR_index = int( TR - y.begin( ) );
            TL_index = int( TM_index - ( TR_index - TM_index ) );
        }
        if( TR_index > int( y.size( ) ) )
        {
            TR_index = int( y.size( ) ) - 1;
        }
        if( TL_index < 0 )
        {
            TL_index = 0;
        }
        peakParams.back( )[1] = x[TM_index];
        peakParams.back( )[3] = TL_index;
        peakParams.back( )[4] = TM_index;
        peakParams.back( )[5] = TR_index;
    }
    
    // loop through peaks, applies activation formula
    for( int i = 0 ; i < int( peakParams.size( ) ) ; i++ )
    {
        peakParams[i][0] = activation( x[peakParams[i][3]], x[peakParams[i][5]],
                                       x[peakParams[i][4]] );
    }
}

// activation formula used in pointsparam
double activation( double TL, double TR, double TM )
{
    double m_g = 0.0, E = 0.0, C = 0.0, K = .000086173303;
    double b = 0.0, t = 0.0, d = 0.0, w = 0.0;
    TL += 273.15;
    TM += 273.15;
    TR += 273.15;
    
    t = TM - TL;
    if( t == 0 )
    {
        t = 1;
    }

    w = TR - TL;
    d = TR - TM;
    m_g = d / w;
    b = 1.58 + 4.2 * ( m_g - 0.42 );
    C = 1.51 + ( 3 * ( m_g - 0.42 ) );

    E = ( ( C * K * ( TM * TM ) ) / t ) - ( b * ( 2 * K * TM ) );
    if( E > 3.0 )
    {
        E = 3;
    }
    return E;
}

void printFindings( std::vector<double>& x, std::vector<double>& y,
                    std::vector<int>& minimum, std::vector<int>& maxima,
                    std::vector<int>& inflectPnt, std::string dir )
{
    std::ofstream myfile;
    myfile.open(dir+"/maxima.csv");
    if(!myfile.is_open()){
        exit(1);
    }
    for(auto i = maxima.begin(); i != maxima.end();++i){
        myfile<<x[*i]<<","<<y[*i]<<"\n";
    }
    myfile.close();
    myfile.open(dir+"/minimum.csv");
    if(!myfile.is_open()){
        exit(1);
    }
    for(auto i = minimum.begin(); i != minimum.end();++i){
        myfile<<x[*i]<<","<<y[*i]<<"\n";
    }
    myfile.close();
    myfile.open(dir+"/inflection.csv");
    if(!myfile.is_open()){
        exit(1);
    }
    for(auto i = inflectPnt.begin(); i != inflectPnt.end();++i){
        myfile<<x[*i]<<","<<y[*i]<<"\n";
    }
    myfile.close();
    myfile.open(dir+"/curve.csv");
    if(!myfile.is_open()){
        exit(1);
    }
    for(int i = 0; i < int(x.size());++i){
        myfile<<x[i]<<","<<y[i]<<"\n";
    }
    myfile.close();
    
}
void write( std::vector<std::vector<double>> glow_curves,
            std::vector<double> y,std::vector<double> x,
            std::string output_name )
{
    std::ofstream file;
    output_name += "/output.csv";
    file.open(output_name);
    if(!file.is_open()){
        exit(1);
    }
    file<<"temp,";
    for(int j = 0; j<int(glow_curves.size());++j){
        std::string ster = "count_" + std::to_string(j);
        file<<","<<ster;
    }
    file<<",\n";
    file.setf(std::ios_base::fixed);
    //file<<std::setprecision(5);
    for(int i = 0; i<int(y.size());++i){
        file << x[i]<<",";
        //file << y[i];
        for(int j = 0; j<int(glow_curves.size());++j){
            file<<","<<double(glow_curves[j][i]);
        }
        file<<",\n";
    }
    file.close();
}

// Curve Fitting for Peaks after Subtraction
void nonMaxPeaks( std::vector<double>& x, std::vector<double>& y,
                  std::vector<double> secDerivative,
                  std::vector<int>& maxima, std::vector<int>& minima,
                  std::vector<std::vector<double>>& peakParams,
                  std::string output_dir)
{
    std::vector<double> yTemp = y;
    std::string dir = output_dir;
    
    // "line sum" of whole curve
    const double origPeakArea = std::accumulate( yTemp.begin( ), yTemp.end( ), 0.0 );
    std::vector<double> sum( x.size( ), 0.0 );
    std::vector<std::vector<double>> peaks;
    // Individual Riemann "Bars" (Integration) for entire curve
    for( int i = 0 ; i < int( peakParams.size( ) ) ; i++ )
    {
        std::vector<double> peak( x.size( ), 0.0 );
        //call FOKModel from FOKMOdel.cpp and populate peak data
        FOKModel( x, peak, peakParams[i][1], peakParams[i][2], peakParams[i][0] );
        //apply plus to all data in peak and store in sum
        transform( peak.begin(), peak.end(), sum.begin(), sum.begin(), std::plus<double>() );
        //substract data in sum sequantially from all sequential yTemp data 
        transform( yTemp.begin(), yTemp.end(), sum.begin(), yTemp.begin(), std::minus<double>() );
    }
    for( int j = 0 ; j < int( yTemp.size( ) ) ; j++ )
    {
        if( yTemp[j] < 0.0 || j > maxima.back( ) )
        {
            yTemp[j] = 0.0;
        }
    }
    std::vector<double> smoothed;
    
    // add the Riemann "Bars" (Integration) over entire curve
    double curPeakArea = std::accumulate( yTemp.begin( ), yTemp.end( ), 0.0 );
    while( curPeakArea > ( 0.2 * origPeakArea ) )
    {
        std::vector<int> remainInflects;
        std::vector<double>::iterator TM = y.begin( );
        if( remainInflects.empty( ) )
        {
            int max = 0;
            double maxVal = 0.0;
            for( int i = 0 ; i < int( yTemp.size( ) ) ; i++ )
            {
                if( yTemp[i] > maxVal )
                {
                    maxVal = yTemp[i];
                    max = i;
                }
            }
            
            // auto tempTM = std::max(yTemp.begin(),yTemp.end());
            TM = y.begin( ) + max;
        }
        else
        {
            double min = 10;
            int minIndex = 0;
            for ( int i = 0 ; i < int( remainInflects.size( ) ) ; i++ )
            {
                if( abs( secDerivative[remainInflects[i]] ) < min )
                {
                    min = abs( secDerivative[remainInflects[i]] );
                    minIndex = remainInflects[i];
                }
            }
            TM = y.begin( ) + minIndex;
        }
        int minIndex = int( TM - y.begin( ) );
        int minLeft = 0, minRight = int( x.size( ) );
        // Look for left adjacent min
        for(int k = minIndex ; k > 0 ; k-- )
        {
            for( auto m = minima.begin( ) ; m != minima.end( ) ; m++ )
            {
                if( k == *m )
                {
                    minLeft = k;
                    break;
                }
            }
            if ( minLeft != 0 )
            {
                break;
            }
        }
        // Look for right adjacent min
        for( int k = minIndex ; k < int( x.size( ) ); k++ )
        {
            for( auto m = minima.begin( ) ; m != minima.end( ) ; m++ )
            {
                if( k == *m )
                {
                    minRight = k;
                    break;
                }
            }
            if ( minRight != int( x.size( ) ) )
            {
                break;
            }
        }
        peakParams.push_back( std::vector<double>( 6, 0.0 ) );
        auto TL = y.begin( );
        auto TR = y.end( );
        // Finds half max points
        double half_intensity = *TM / 2.0;
        TL = lower_bound( TL + minLeft, TM, half_intensity, std::less<double>( ) );
        //TR = lower_bound( TM, TR + minRight, half_intensity, std::greater<double>( ) );
        TR = lower_bound(TM, TR, half_intensity, std::greater<double>());
        // Calculates half width half max
        int diff1 = int( TM - TL );
        int diff2 = int( TR - TM );
        if( TR == y.end( ) )
        {
            diff2 += diff1;
        }
        int TM_index = int( TM - y.begin( ) );
        int TL_index, TR_index;
        // Centering, finds smaller out of left and right side and adjust accordingly
        if( diff1 <= diff2 || *TR > *TM )
        {
            TL_index = int( TL - y.begin( ) );
            TR_index = int( TM_index + ( TM_index - TL_index ) );
        }
        else if( *TL > *TM )
        {
            TR_index = int( TR - y.begin( ) );
            TL_index = int( TM_index - ( TR_index - TM_index ) );
        }
        else
        {
            TR_index = int( TR - y.begin( ) );
            TL_index = int( TM_index - ( TR_index - TM_index ) );
        }
        if( TR_index > int( y.size( ) ) )
        {
            TR_index = int( y.size( ) ) - 1;
        }
        if( TL_index < 0 )
        {
            TL_index = 0;
        }
        peakParams.back( )[1] = x[TM_index];
        peakParams.back( )[2] = y[TM_index];
        peakParams.back( )[3] = TL_index;
        peakParams.back( )[4] = TM_index;
        peakParams.back( )[5] = TR_index;
        peakParams.back( )[0] = activation( x[TL_index], x[TR_index], x[TM_index] );
        
        std::vector<double> peak( x.size( ), 0.0 );
        FOKModel( x, peak, x[TM_index], y[TM_index], peakParams.back( )[0] );
        transform( peak.begin( ), peak.end( ), sum.begin() , sum.begin( ),
                   std::plus<double>( ) );
        transform( yTemp.begin( ), yTemp.end( ), sum.begin( ), yTemp.begin( ),
                   std::minus<double>( ) );
        curPeakArea = std::accumulate( yTemp.begin( ), yTemp.end( ), 0.0 );
    }
}
