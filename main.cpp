#include <windows.h>

#include <iostream>
#include <vector>

#include "matrix.h"
using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::vector;

vector<double> multiply_matrix_by_vector(Matrix M, vector<double> vec)
{
    vector<double> multiplied_vector(vec.size());
    for (int i = 0; i < vec.size(); ++i)
    {
        double sum = 0;
        for (int j = 0; j < vec.size(); ++j)
        {
            sum += M.value[i][j] * vec[i];
        }

        multiplied_vector[i] = sum;
    }

    return multiplied_vector;
}

// this program on top of finding the inverse matrix, solves system Ax=b. If you
// are not interested in x vector, just look at the inverse matrix output, since further
// steps don't affect inverse matrix calculation
int main()
{
    // file matrix.in contains the number of rows, vector b (in case you want to solve system Ax=b)
    freopen("matrix.in", "r", stdin);
    // uncomment the row below to get the output in a file
    // freopen("inverse.out", "w", stdout);
    int size;
    cin >> size;
    vector<double> b(size);
    for (int i = 0; i < size; ++i)
        cin >> b[i];

    Matrix A(size, size);
    for (int i = 0; i < A.get_rows(); ++i)
        for (int j = 0; j < A.get_cols(); ++j)
            cin >> A.value[i][j];

    if (A.det() != 0)
    {
        A.add_identical_matrix();
        A.gauss_jordan();
        A.extract_inverse_matrix();
        cout << "inverse matrix: \n";
        A.print();
        vector<double> x_is = multiply_matrix_by_vector(A, b);

        cout << "matrix is invertable; there is one exact solution:x = A^(-1)*b" << endl
             << "x: ";

        for (auto &element : x_is)
            cout << element << '\t';
    }
    else
    {
        A.extend_matrix_with_vector(b);
        A.gauss_jordan();
        cout << "matrix A is singular, can't find the inverse" << endl;

        if (A.value[A.get_rows() - 1][A.get_cols() - 1] != 0)
        {
            cout << "no solutions\n";
        }
        else
        {
            cout << "infenitely many solutions\n";

            // A.print_vec_inf_solutions();
        }
    }
}
