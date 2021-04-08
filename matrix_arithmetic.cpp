//
//  matrix_arithmetic.cpp
//  GlowCurveAnalsys
//
//  Initially created by jeremy hepker on 3/31/19.
//
//  Modified and re-organized by Jack Yu UROP 2020 Fall
//

#include "Levenberg_Marquardt.hpp"

using namespace std;
//---------------------Levenburg Marquardt Method Helper Functions----------------------//

double First_Order_Kinetics::Deriv2(const double input, const vector<double> params, int n)
{
    // Assumes input is a single row matrix
    double DERIV_STEP = (params[n] < 5) ? 1e-2:0.01;
    // Returns the derivative of the nth parameter
    vector<double> params1 = params;
    vector<double> params2 = params;
    double d = 0;
    // Use central difference to get derivative
    params1[n] -= DERIV_STEP;
    params2[n] += DERIV_STEP;
    double p1 = Func2(input, params1);
    double p2 = Func2(input, params2);
    d = (p2 - p1) / (2*DERIV_STEP);
    return d;
}

//Function to transpose a matrix
void First_Order_Kinetics::transpose(vector<vector<double>> const &A,vector<vector<double>> &B, int n, int m){
    for(auto i = B.begin(); i != B.end();++i) i->resize(A.size());
    for (int i=0; i < n; i++) {
        for (int j=0; j < m; j++) {
            B[j][i] = A[i][j];
        }
    }
}
//Functtion to multiply two matricies together
vector<vector<double>> First_Order_Kinetics::multiply(vector<vector<double>> const &A,vector<vector<double>> const &B){
    size_t row_A = A.size();
    size_t col_A = A[0].size();
    size_t col_B = B[0].size();
    double stop = 1* pow(double(10),double(-35));
    vector<vector<double>> C(row_A, vector<double>(col_B,0));
    for (size_t i = 0; i <row_A; i++){
        for (size_t j = 0; j < col_B; j++){
            for (size_t k = 0; k < col_A; k++){
                C[i][j] += A[i][k] * B[k][j];
                ((C[i][j] < stop) ? C[i][j] = 0.0 :false);
            }
        }
    }
    return C;
};
//Function to muliply an array with a matrix
vector<double> First_Order_Kinetics::vec_matrix_multi(vector<vector<double>> const &A,vector<double> const &B){
    vector<double> output;
    for(auto i = A.begin(); i != A.end();++i){
        double sum = 0.0;
        for(int j = 0; j != int(B.size()); ++j){
            sum += i->at(j) * B[j];
        }
        output.push_back(sum);
    }
    return output;
}
//Function to invert a matrix, with option for negtive inverse.
void First_Order_Kinetics::invert(vector<vector<double>> &A, bool neg){
    if(A.size() == 2){
        double det = A[0][0] *A[1][1] - A[1][0] * A[0][1];
        vector<vector<double>> temp(2, vector<double>(2,0));
        if(neg){
            det *= (-1.0);
        }
        temp[0][0] = (1/det)*A[1][1];
        temp[0][1] = -(1/det)*A[0][1];
        temp[1][0] = -(1/det)*A[1][0];
        temp[1][1] = (1/det)*A[0][0];
        A = temp;
    }else{
        int N = int(A.size());
        double det = determinant(A, N);
        if (det == 0)
        {
            throw det;
        }
        
        // Find adjoint
        vector<vector<double>> adj(A.size(), vector<double>(A.size(), 0.0));
        adjoint(A, adj);
        
        // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
        for (int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                A[i][j] = adj[i][j]/det;
            }
        }
        
        
    }
}
void First_Order_Kinetics::adjoint(vector<vector<double>> &A,vector<vector<double>> &adj){
    int N = int(A.size());
    if (N == 1){
        adj[0][0] = 1;
        return;
    }
    
    // temp is used to store cofactors of A[][]
    double sign = 1.0;
    vector<vector<double>> temp(N, vector<double>(N,0.0));
    
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            // Get cofactor of A[i][j]
            cofactor(A, temp, i, j, N);
            
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1.0: -1.0;
            
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign)*(determinant(temp, N-1));
        }
    }
}
void First_Order_Kinetics::cofactor(vector<vector<double>> &A, vector<vector<double>> &temp, int p, int q, int n){
    int i = 0, j = 0;
    
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++){
        for (int col = 0; col < n; col++){
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q){
                temp[i][j++] = A[row][col];
                
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1){
                    j = 0;
                    i++;
                }
            }
        }
    }
}

double First_Order_Kinetics::determinant(vector<vector<double>> &A, int size){
    double s=1.0,det=0.0;
    vector<vector<double>> m_minor(size, vector<double>(size,0.0));
    double i,j,m,n,c;
    if (size==1){
        return (A[0][0]);
    }
    else{
        det=0;
        for (c=0;c<size;c++){
            m=0;
            n=0;
            for (i=0;i<size;i++){
                for (j=0;j<size;j++){
                    m_minor[i][j]=0;
                    if (i != 0 && j != c){
                        m_minor[m][n]=A[i][j];
                        if (n<(size-2)) n++;
                        else{
                            n=0;
                            m++;
                        }
                    }
                }
            }
            det=det + s * (A[0][c] * determinant(m_minor,size-1));
            s=-1.0 * s;
        }
    }
    return (det);
}
//Functin to preform the dot product between two vectors
double First_Order_Kinetics::dotProduct(vector<double> A, vector<double> B){
    double product = 0;
    // Loop for calculate cot product
    for (int i = 0; i < int(A.size()); i++)
    product = product + A[i] * B[i];
    return product;
}
//Function to create an identity matrix multiplied by give lambda
vector<vector<double>> First_Order_Kinetics::Identity(int num, double lambda){
    int row, col;
    vector<vector<double>> I(num, vector<double>(num,0));
    for (row = 0; row < num; row++){
        for (col = 0; col < num; col++){
            // Checking if row is equal to column
            if (row == col)
            I[row][col] = lambda;
        }
    }
    return I;
}

// Calculate First Derivatives with temperature and data input
void First_Order_Kinetics::deriv(std::vector<double>& x, std::vector<double>& y, std::vector<double>& derivative)
{
    derivative.resize(x.size());
    int size = int(x.size()) - 2;
    // Five Point Stencil Method
    for (int i = 2; i < size; i++)
    {
        derivative[i - 2] = (-y[i + 2] + 8.0 * y[i + 1] - 8.0 * y[i - 1] + y[i - 2]) / 12.0;
    }
    size = int(derivative.size()) - 1;
    int sign = 4, lastSign = 0;
    bool positive = true;
    for (int i = 0; i < size; i++)
    {
        if (derivative[i] < 0.0 && positive && sign >= 4)
        {
            positive = false;
            lastSign = sign;
            sign = 1;
        }
        else if (derivative[i] < 0.0 && positive && sign < 4)
        {
            for (int j = 1; j <= sign; j++)
            {
                derivative[i - j] = -(derivative[i - j]);
                positive = false;
            }
            sign = lastSign + sign;
        }
        else if (derivative[i] > 0.0 && !positive && sign < 4)
        {
            for (int j = 1; j <= sign; j++)
            {
                derivative[i - j] = abs(derivative[i - j]);
                positive = true;
            }
            sign = lastSign + sign;
        }
        else if (derivative[i] > 0.0 && !positive && sign >= 4)
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
