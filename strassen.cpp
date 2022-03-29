#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <cstdarg>

using namespace std;

//This crossover point was found theoretically
const int CROSS_OVER_POINT = 15;

//Initialize random number generator for matrix values between 0 and 1
random_device dev;
mt19937 mersenne(dev());
uniform_int_distribution<mt19937::result_type> value_gen(0,2);

// Function to initialize matrices that will be multiplies


vector<int> matrix_add(vector<int> matrix_1, vector<int> matrix_2){
    
    vector<int> sum_matrix;

    /* Iterate through both input matrices and store the sum
    of corresponding elements in the corresponding position in
    sum matrix */
    for (int i = 0; i < matrix_1.size(); i++){
        sum_matrix.push_back(matrix_1[i] + matrix_2[i]);
    }

    return sum_matrix;
}


vector<int> matrix_sub(vector<int> matrix_1, vector<int> matrix_2){

    vector<int> diff_matrix;

    for (int i = 0; i < matrix_1.size(); i++){
        diff_matrix.push_back(matrix_1[i] - matrix_2[i]);
    }

    return diff_matrix;
}

vector<vector<int>> matrix_splitter(vector<int> matrix, int n){

    //Divide the input matrices into 4 parts, one iteration for each part
    //and place all 4 parts into a super list

    vector<int> quad_1;
    vector<int> quad_2;
    vector<int> quad_3;
    vector<int> quad_4;

    int offset = n/2;

    //Create first quadrant
    for (int row = 0; row < n/2; row++){
        for (int col = 0; col < n/2; col++){
            quad_1.push_back(matrix[col+(row*n)]);
        }
    }

    //Create second quadrant
    for (int row = 0; row < n/2; row++){
        for (int col = 0; col < n/2; col++){
            quad_2.push_back(matrix[(col+offset)+(row*n)]);
        }
    }

    //Create third quadrant
    for (int row = 0; row < n/2; row++){
        for (int col = 0; col < n/2; col++){
            quad_3.push_back(matrix[col+((row + offset)*n)]);
        }
    }

    //Create fourth quadrant 
    for (int row = 0; row < n/2; row++){
        for (int col = 0; col < n/2; col++){
            quad_4.push_back(matrix[(col+offset)+((row + offset)*n)]);
        } 
    }

    vector<vector<int>> quad_list = {quad_1 ,quad_2, quad_3, quad_4};

    return quad_list;

}

vector<int> standard_mm(vector<int> factor_1, vector<int> factor_2, int n){

    vector<int> product;
    int entry_sum = 0;

    // Iterates through the n^2 positions in the product matrix
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++){
            /* factor_1 indices will traverse k entries on a row for each j
                factor_2 indices will traverse k entries on a column for each j,
                 */
                entry_sum += factor_1[k + n*i] * factor_2[n*k + j];
            }
            product.push_back(entry_sum);
            entry_sum = 0;
        }
        
    }

    return product;

}

vector<int> strassens_mm(vector<int> factor_1, vector<int> factor_2, int n){

    vector<int> product;
    //Switch to regular matrix multiplication when we cross threshold
    if (n <= CROSS_OVER_POINT){
        return standard_mm(factor_1, factor_2, n);
    }

    vector<vector<int>> quads_1 = matrix_splitter(factor_1, n);
    vector<vector<int>> quads_2 = matrix_splitter(factor_2, n);
    
    //Make 7 recursive calls for P1-P7
    vector<int> P1 = strassens_mm(quads_1[0], matrix_sub(quads_2[1], quads_2[3]), n/2);
    vector<int> P2 = strassens_mm(matrix_add(quads_1[0], quads_1[1]), quads_2[3], n/2);
    vector<int> P3 = strassens_mm(matrix_add(quads_1[2], quads_1[3]), quads_2[0], n/2);
    vector<int> P4 = strassens_mm(quads_1[3], matrix_sub(quads_2[2], quads_2[0]), n/2);
    vector<int> P5 = strassens_mm(matrix_add(quads_1[0], quads_1[3]), matrix_add(quads_2[0], quads_2[3]), n/2);
    vector<int> P6 = strassens_mm(matrix_sub(quads_1[1], quads_1[3]), matrix_add(quads_2[2], quads_2[3]), n/2);
    vector<int> P7 = strassens_mm(matrix_sub(quads_1[2], quads_1[0]), matrix_add(quads_2[0], quads_2[1]), n/2);

    /*Combine the values returned from P1-P7 into the 
    4 terms needed in the prdouct matrix */

    vector<int> term_A = matrix_add(matrix_add(P5,P6), matrix_sub(P4,P2));
    vector<int> term_B = matrix_add(P1,P2);
    vector<int> term_C = matrix_add(P3,P4);
    vector<int> term_D = matrix_add(matrix_sub(P1, P3), matrix_add(P5,P7));
    
    //**Reconstruct the product matrix from the four terms***
    for (int row = 0; row < n; row++){

        //Build top half of complete product matrix
        if (row < n/2){
            for (int col = 0; col < n/2; col++){
                product.push_back(term_A[col+(row*(n/2))]);
            }
            for (int col = 0; col < n/2; col++){
                product.push_back(term_B[col+(row*(n/2))]);
            }

        //Build bottom half of product matrix with last 2 terms
        } else {
            for (int col = 0; col < n/2; col++){
                product.push_back(term_C[col+(row-(n/2))*(n/2)]);
            }
            for (int col = 0; col < n/2; col++){
                product.push_back(term_D[col+(row-(n/2))*(n/2)]);   
            }
        }
    }
    return product;
}

int main(int argc, char** argv){

    vector<int> product;
    string line;

    int n = atoi(argv[2]);
    int next_pow2 = pow(2, ceil(log2(n)));

    vector<int> factor_1 = matrix_constructor(n);
    vector<int> factor_2 = matrix_constructor(n);
    int file_itr = 0;
    int matrix_itr = 0;

    while(getline(input_file, line)){
        if (file_itr < pow(n,2)){
            //Add padding to input matrix at the end of every row as necessary
            factor_1.push_back(line);
            mat_itr++;
            file_itr++;

            if ((file_itr % (n-1) == 0) && n != next_pow2){
                for (int i = 0; i < n-next_pow2; i++){
                    factor_1.push_back(0);
                    mat_itr++;
                }
            }

            //Pad the remaining positions in the matrix until its dimensions are the nearest power of 2
            if (file_itr == pow(n,2) - 1 && pow(n,2) != pow(next_pow2, 2)){
                while mat_itr < pow(next_pow2, 2)){
                    factor_1.push_back(0);
                    mat_itr++;
                }
                mat_itr = 0
            }

        } else {
            factor_2.push_back(line);
            mat_itr++;
            file_itr++;

            if ((file_itr % (n-1) == 0) && n != next_pow2){
                for (int i = 0; i < n-next_pow2; i++){
                    factor_2.push_back(0);
                    mat_itr++;
                }
            }

            //Pad the remaining positions in the matrix until its dimensions are the nearest power of 2
            if (file_itr == (2*pow(n,2)) - 1 && pow(n,2) != pow(next_pow2, 2)){
                while mat_itr < pow(next_pow2, 2)){
                    factor_2.push_back(0);
                    mat_itr++
                }
            }
    
    }

    for (int i = 0; i < factor_1.size(); i++){
        cout << factor_1[i] << " ";
    }
    cout << endl;

    for (int i = 0; i < factor_2.size(); i++){
        cout << factor_2[i] << " ";
    }
    cout << endl;

    if (n <= CROSS_OVER_POINT){
        product = standard_mm(factor_1, factor_2, sqrt(factor_1.size()));
    } else { 
        product = strassens_mm(factor_1, factor_2, sqrt(factor_1.size()));
    }

    for (int i = 0; i < next_pow2; i++){
        cout << product[i + (i*next_pow2)];
        cout << " ";
    }   
}


vector<int> matrix_constructor(int n){

    vector<int> matrix;
    int next_pow2 = pow(2,ceil(log2(n)));

    //Different rows will be indexed by offset of n*(row_num-1)
        for (int i = 0; i < next_pow2; i++){

            //Fill in the first n columns of the first n rows with random entries
            if (i < n){
                for (int nums = 0; nums < n; nums++){
                    matrix.push_back(value_gen(mersenne));
                }
                //Fill in every column after the nth with padding (0's)
                for (int pad = 0; pad < (next_pow2 - n); pad++){
                    matrix.push_back(0);
                }

            //Fill in any row below the nth row with padding (0's)
            } else {
                for (int remn = 0; remn < next_pow2; remn++){
                    matrix.push_back(0);
                }
            }   
        }
    return matrix;

}