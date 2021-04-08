//
//	Savitzsky_Golay.hpp
//  GlowCurveAnalsys
//
//  Created by Jack YuUROP 2020 Fall
//  Matrix aithmetics written by Jeremy Hepker
//
#include "Savitzky_Golay.hpp"

void SG_smooth(vector<double>& y, int window_size, int size) {
    //create the Vandermonde Matrix
    vector<vector<double>> A(window_size, vector<double>(size));
    int up = -(window_size / 2);
    int down = 1;
    //denote the degree of power
    int rate = 0;
    //fill in the upper half of the matrix
    for (int i = 0; i < window_size / 2; i++) {
        for (int j = 0; j < size; j++) {
            A[i][j] = pow(up, rate);
            rate++;
        }
        up++;
        rate = 0;
    }
    //fill in the middle vector with leading 1
    A[window_size / 2][0] = 1;
    for (int k = 1; k < size; k++)
        A[window_size / 2][k] = 0;
    //fill in the lower half of the matrix
    for (int l = window_size / 2 + 1; l < window_size; l++) {
        for (int m = 0; m < size; m++) {
            A[l][m] = pow(down, rate);
            rate++;
        }
        down++;
        rate = 0;
    }
    //compute the transpose of A
    vector<vector<double>> A_T(size, vector<double>(window_size));
    transpose(A, A_T, window_size, size);
    //compute the inverse of the multiplication of A_T and A
    vector<vector<double>> mul = multiply(A_T, A);
    invert(mul, false);
    //multiply the inverse with the tranpose to get the multiplication matrix
    vector<vector<double>> M = multiply(mul, A_T);
    //take the first row of M and save as change
    vector<double> change = M[0];
    
    //extend the left and right side of count data with the reverse of the original data
    vector<double> y_temp(y.size() + window_size - 1);
    for (int i = 0; i < window_size / 2; i++) {
        y_temp[i] = y[window_size / 2 - i - 1];
    }
    for (int i = 0; i < static_cast<int>(y.size()); i++) {
        y_temp[i + window_size / 2] = y[i];
    }
    int index = 1;
    for (int i = static_cast<int>(y.size()) + window_size / 2; i < static_cast<int>(y_temp.size()); i++) {
        y_temp[i] = y_temp[i - index];
        index += 2;
    }
    
    //convolve y_temp with the change matrix 
    vector<double> y_new(y_temp.size() - window_size + 1);
    int i = window_size - 3;
    while (i < static_cast<int>(y_temp.size()) - 1) {
        for (int j = 0; j < static_cast<int>(change.size()) - 1; j++) {
            if (i > window_size - 3) {
                y_new[i - window_size + 2] += change[j] * y_temp[i - j];
            }
        }
        i += 1;
    }
    
    y = y_new;

    //for (int i = 0; i < 9; i++) {
    //    for (int j = 0; j < 9; j++) {
    //        cout << m[i][j] << " ";
    //    }
    //    cout << endl;
    //}
}

//take the transpose of a matrix
void transpose(vector<vector<double>> const& A, vector<vector<double>>& B, int n, int m) {
    //for (auto i = B.begin(); i != B.end(); ++i) i->resize(A.size());
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            B[j][i] = A[i][j];
        }
    }
}
//Functtion to multiply two matricies together
vector<vector<double>> multiply(vector<vector<double>> const& A, vector<vector<double>> const& B) {
    size_t row_A = A.size();
    size_t col_A = A[0].size();
    size_t col_B = B[0].size();
    vector<vector<double>> C(row_A, vector<double>(col_B, 0));
    for (size_t i = 0; i < row_A; i++) {
        for (size_t j = 0; j < col_B; j++) {
            for (size_t k = 0; k < col_A; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

double determinant(vector<vector<double>>& A, int size) {
    double s = 1.0, det = 0.0;
    vector<vector<double>> m_minor(size, vector<double>(size, 0.0));
    double i, j, m, n, c;
    if (size == 1) {
        return (A[0][0]);
    }
    else {
        det = 0;
        for (c = 0; c < size; c++) {
            m = 0;
            n = 0;
            for (i = 0; i < size; i++) {
                for (j = 0; j < size; j++) {
                    m_minor[i][j] = 0;
                    if (i != 0 && j != c) {
                        m_minor[m][n] = A[i][j];
                        if (n < (size - 2)) n++;
                        else {
                            n = 0;
                            m++;
                        }
                    }
                }
            }
            det = det + s * (A[0][c] * determinant(m_minor, size - 1));
            s = -1.0 * s;
        }
    }
    return (det);
    //cout << "finish determinant" << endl;
}

void cofactor(vector<vector<double>>& A, vector<vector<double>>& temp, int p, int q, int n) {
    int i = 0, j = 0;

    // Looping for each element of the matrix
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q) {
                temp[i][j++] = A[row][col];
                //for (int a = 0; a < 6; a++) {
                //    for (int b = 0; b < 6; b++) {
                //        cout << temp[a][b] << " ";
                //    }
                //    cout << endl;
                //}
                //cout << endl;
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
    //cout << "finish cofactor" << endl;
}

void adjoint(vector<vector<double>>& A, vector<vector<double>>& adj) {
    int N = int(A.size());
    if (N == 1) {
        adj[0][0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    double sign = 1.0;
    vector<vector<double>> temp(N, vector<double>(N, 0.0));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // Get cofactor of A[i][j]
            cofactor(A, temp, i, j, N);
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(temp, N - 1));
            //for (int a = 0; a < 6; a++) {
            //    for (int b = 0; b < 6; b++) {
            //        cout << adj[a][b] << " ";
            //    }
            //    cout << endl;
            //}
            //cout << endl;
        }
    }
    //cout << "finish adjoint" << endl;
}

//Function to invert a matrix, with option for negtive inverse.
void invert(vector<vector<double>>& A, bool neg) {
    if (A.size() == 2) {
        double det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
        vector<vector<double>> temp(2, vector<double>(2, 0));
        if (neg) {
            det *= (-1.0);
        }
        temp[0][0] = (1 / det) * A[1][1];
        temp[0][1] = -(1 / det) * A[0][1];
        temp[1][0] = -(1 / det) * A[1][0];
        temp[1][1] = (1 / det) * A[0][0];
        A = temp;
    }
    else {
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
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = adj[i][j] / det;
            }
        }
    }
    //cout << "finish invert" << endl;
}

