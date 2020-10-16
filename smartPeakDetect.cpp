//
//  smartPeakDetect.cpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 7/9/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#include "smartPeakDetect.hpp"

// Curve Fitting and Vetting
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
    
    // Pushes all current Maxes and Mins into vector
    for( int i = 1 ; i < size ; i++ )
    {
        if( derivative[i] < 0.0 && derivative[i - 1] > 0.0 )
        {
            maxima.push_back( i );
        }
        if( derivative[i] > 0.0 && derivative[i - 1] < 0.0 )
        {
            minimum.push_back( i );
        }
    }
    minimum.push_back( int( x.size( ) ) - 1 );
    
    size = int( secDerivative.size( ) );
    // Pushes inflection points into vector
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
    
    // *************** //
    // TODO::
    // *************** //
    
    // Original Value = 10
    int adj_peak_distance = 10;
    // Original Value = 5
    int near_point_height = 5;
    // Original Value = 10
    // int low_threshold = 0;
    
    // loop through maximas
    for( int i = 1 ; i < int( maxima.size( ) ) - 1 ; i++ )
    {
        // if peak is not as tall as adjacent peaks
        if( y[maxima[i]] < y[maxima[i + 1]] && y[maxima[i]] < y[maxima[i - 1]] )
        {
            maxima.erase( maxima.begin( ) + i );
            i--;
            continue;
        }
        
        // if adjacent peaks are not far enough from each other
        if( abs( x[maxima[i]] - x[maxima[i - 1]] ) <= adj_peak_distance )
        {
            if( abs( x[maxima[i]] - x[maxima[i + 1]] ) <= adj_peak_distance )
            {
                maxima.erase( maxima.begin( ) + i + 1 );
                maxima.erase( maxima.begin( ) + i - 1 );
                i--;
                continue;
            }
            else
            {
                maxima.erase( maxima.begin( ) + i );
                i--;
                continue;
            } // else
        } // if
    } // for
    
    /*
    // loop through all peaks
    for( int i = 0 ; i < int( maxima.size( ) ) ; i++ )
    {
        // if a peak is not as tall as a previous peak or is less than a certain threshold
        if( ( y[maxima[i]] < y[maxima[i] - near_point_height] ) ||
              y[maxima[i]] < low_threshold )
        {
            maxima.erase( maxima.begin( ) + i );
            i--;
            continue;
        }
    }
    */
    
    // loop through all peaks
    for( int i = 0 ; i < int( maxima.size( ) ) ; i++ )
    {
        // if a peak is not as tall as a previous peak or is less than a certain threshold
        if( ( y[maxima[i]] < y[maxima[i] - near_point_height] ) )
        {
            maxima.erase( maxima.begin( ) + i );
            i--;
            continue;
        }
    }
    
    // *************** //
    // *************** //
    // *************** //
}

// Full Width Half Max
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
    
    // loop through peaks, applies formula
    for( int i = 0 ; i < int( peakParams.size( ) ) ; i++ )
    {
        peakParams[i][0] = activation( x[peakParams[i][3]], x[peakParams[i][5]],
                                       x[peakParams[i][4]] );
    }
}

// Formula
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

// Sudo Main function for smartPeakDetect
void findPeaks( std::vector<double>& x, std::vector<double>& y,
                std::vector<std::vector<double>>& peakParams )
{
    std::vector<double> xNew = x, yNew= y;
    std::vector<std::vector<double>> peaks;
    std::vector<int> maximums, minimum, inflections;
    std::vector<double> firstDir( x.size( ), 0.0 );
    std::vector<double> secDir( x.size( ), 0.0 );
    firstDeriv( xNew, yNew, firstDir );
    secDeriv( xNew, yNew, secDir );
    smartPoints( xNew, yNew, minimum, maximums,firstDir, secDir, inflections );
    pointsParams( xNew, yNew, maximums, minimum, peakParams );
    // Copy result into peakParams
    for( int i = 0 ; i < int( peakParams.size( ) ) ; i++ )
    {
        peakParams[i][2] = yNew[peakParams[i][4]];
    }
    nonMaxPeaks( xNew, yNew, secDir,maximums, minimum, peakParams );
    
    /*THIS IS USED FOR OUTPUTING PEAK DETECTION STATS JUST UNCOMMENT*/
    minimum.clear( );
    for( int j = 0 ; j < int( peakParams.size( ) ) ; j++ )
    {
        minimum.push_back( peakParams[j][3] );
        minimum.push_back( peakParams[j][5] );
    }
    //change to local path need fix!
    printFindings( xNew, yNew, minimum, maximums, inflections,
                "/Users/rhellab/Desktop/GCA/GlowCurveAnalysis/" );
    for( auto i = peakParams.begin( ) ; i != peakParams.end( ) ; i++ )
    {
        std::vector<double> peak( x.size( ), 0.0 );
        FOKModel( xNew, peak, i->at( 1 ), i->at( 2 ), i->at( 0 ) );
        peaks.push_back( peak );
    }
    write( peaks, yNew, xNew, "/Users/rhellab/Desktop/GCA/GlowCurveAnalysis/" );
}

// Calculate First Derivatives
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

// Calculate Second Derivatives
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
    output_name += "_output.csv";
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
                  std::vector<std::vector<double>>& peakParams )
{
    std::vector<double> yTemp = y;
    //changed to local path
    std::string dir = "/Users/rhellab/Desktop/GCA/GlowCurveAnalysis/";
    
    // *************** //
    // int iteration = 0;
    // *************** //
    
    // "line sum" of whole curve
    const double origPeakArea = std::accumulate( yTemp.begin( ), yTemp.end( ), 0.0 );
    std::vector<double> sum( x.size( ), 0.0 );
    std::vector<std::vector<double>> peaks;
    // Individual Riemann "Bars" (Integration) for entire curve
    for( int i = 0 ; i < int( peakParams.size( ) ) ; i++ )
    {
        std::vector<double> peak( x.size( ), 0.0 );
        FOKModel( x, peak,peakParams[i][1], peakParams[i][2], peakParams[i][0] );
        
        /*
        // *************** //
        // *************** //
        // *************** //
        
        std::ofstream out;
        std::stringstream sstream;
        std::string out_name = "/Users/raymondlu/Downloads/GlowCurveAnalsys-master_new/GlowCurveAnalysis/nonMaxPeak_";
        
        sstream << iteration;
        out_name += sstream.str( );
        out_name += "_output.csv";
        
        
        out.open(out_name);
        if(!out.is_open()){
            exit(1);
        }
        
        out << "yTemp,peak";
        out <<",\n";
        for(int j = 0; j<int(x.size());++j){
            out << yTemp[j] << ",";
            out << peak[j];
            out <<",\n";
        }
        out.setf(std::ios_base::fixed);
        out<<std::setprecision(5);
        
        std::cout << "Output File : "<< out_name << "\n";
        
        out.close();
        
        iteration++;
        
        // *************** //
        // *************** //
        // *************** //
        */
        
        transform( peak.begin( ), peak.end( ), sum.begin( ), sum.begin( ),
                   std::plus<double>( ) );
        transform( yTemp.begin( ), yTemp.end( ), sum.begin( ), yTemp.begin( ),
                   std::minus<double>( ) );
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
        TR = lower_bound( TM, TR + minRight, half_intensity, std::greater<double>( ) );
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
