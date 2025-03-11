#include "matrix.h"

#include <windows.h>

#include <iostream>

using std::cin;
using std::cout;
using std::endl;
using std::vector;

double round_two_points(double a)
{
    return double(int(a * 100)) / 100.0;
}
// overloads
template <typename T>
vector<T> operator+(const vector<T> &v1, const vector<T> &v2)
{
    vector<T> result(v1.size());

    for (int i = 0; i < v1.size(); ++i)
        result[i] = v1[i] + v2[i];
    return result;
}

template <typename T>
vector<T> operator-(const vector<T> &v1, const vector<T> &v2)
{
    vector<T> result(v2.size());

    for (int i = 0; i < v1.size(); ++i)
        result[i] = v1[i] + v2[i];
    return result;
}

template <typename T> // vector by vector
vector<T> operator*(const vector<T> &v1, const vector<T> &v2)
{
    vector<T> result(v1.size());

    for (int i = 0; i < v1.size(); ++i)
        result[i] = v1[i] * v2[i];
    return result;
}

template <typename t> // scalar by vector
vector<t> operator*(const vector<t> &v1, const t &scalar)
{
    vector<t> result(v1.size());

    for (int i = 0; i < v1.size(); ++i)
        result[i] = v1[i] * scalar;
    return result;
}

template <typename T> // vector by scalar
vector<T> operator*(const T &scalar, const vector<T> &v1)
{
    vector<T> result(v1.size());

    for (int i = 0; i < v1.size(); ++i)
        result[i] = v1[i] * scalar;
    return result;
}

template <typename T>
vector<T> operator/(const vector<T> &v1, const T &scalar)
{
    vector<T> result(v1.size());

    for (int i = 0; i < v1.size(); ++i)
        result[i] = v1[i] / scalar;
    return result;
}

// end of overloads

int Matrix::get_rows()
{
    return rows;
}

void Matrix::set_rows(int n)
{
    this->rows = n;
}

int Matrix::get_cols()
{
    return cols;
}

void Matrix::set_cols(int n)
{
    this->cols = n;
}

void Matrix::print()
{
    for (int i = 0; i < get_rows(); ++i)
    {
        for (int j = 0; j < get_cols(); ++j)
        {
            cout << value[i][j] << "\t";
        }

        cout << "\n";
    }
}

// returns a matrix without a selected row and a column
vector<vector<double>> Matrix::get_minor(int deleted_row, int deleted_column)
{
    vector<vector<double>> minor_matrix(rows - 1, vector<double>(cols - 1));

    for (int i = 0, minor_i = 0; i < rows; ++i)
    {
        if (i == deleted_row)
            continue;

        for (int j = 0, minor_j = 0; j < cols; ++j)
        {
            if (j == deleted_column)
                continue;

            minor_matrix[minor_i][minor_j] = value[i][j];

            ++minor_j;
        }
        ++minor_i;
    }

    return minor_matrix;
}

double Matrix::det()
{
    int result = 0;

    if (value.size() == 1)
        return value[0][0];
    if (value.size() == 2)
        return value[0][0] * value[1][1] - value[0][1] * value[1][0]; // ad-bc

    int sign = 1;

    for (int i = 0; i < get_cols(); ++i)
    {
        Matrix minor(rows - 1, cols - 1);
        minor.value = get_minor(0, i);

        result += sign * value[0][i] * minor.det();
        sign *= -1;
    }

    return (int)result;
}

void Matrix::swap_rows(int a, int b)
{
    int temporary_element;

    for (int i = 0; i < cols; ++i)
    {
        temporary_element = value[a][i];

        value[a][i] = value[b][i];
        value[b][i] = temporary_element;
    }
}

vector<double> Matrix::add_rows_by_index(int a, int b)
{
    vector<double> temporary_row(cols);

    for (int i = 0; i < cols; i++)
        temporary_row[i] = value[a][i] + value[b][i];

    return temporary_row;
}

// extends matrix A with identical matrix I_n to the right A-> (A|I)
void Matrix::add_identical_matrix()
{
    set_cols(cols + cols);
    for (int i = 0; i < rows; ++i)
        value[i].resize(cols);

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
        {
            if (j < int((cols + 1) / 2))
                continue;

            value[i][j] = int(i == j - rows);
        }
}

void Matrix::replace_row(int row_number, vector<double> replacement_row)
{ // todo
    for (int i = 0; i < cols; ++i)
        value[row_number][i] = round_two_points(replacement_row[i]);
}
// (A|I)->(I|A^(-1))
void Matrix::gauss_jordan()
{
    print();
    if (value[0][0] == 0.0) // in making sure the row is not multiplied be 0 in the next loop
        for (int i = 0; i < rows; ++i)
            if (value[i][0] != 0)
            {
                replace_row(0, value[0] + value[i]);
                break;
            }

    for (int i = 0; i < rows; ++i)
    {
        if (value[i][i] == 0)
        {
            bool swapped = false;
            for (int j = i + 1; j < rows; ++j)
            {
                if (value[j][i] != 0)
                {
                    swap_rows(i, j);
                    swapped = true;
                    break;
                }
            }

            if (!swapped)
                continue;
        };

        for (int j = 0; j < rows; ++j)
        {
            if (i == j)
                continue; // we don't want to multiply the main diagonal

            // making sure we get correct sign, so we get zeros everywhere except the main diagonal
            int sign =
                ((value[j][i] > 0 and value[i][i] > 0) or
                 (value[j][i] < 0 and value[i][i] < 0))
                    ? sign = -1
                    : sign = 1;

            double modified_value = abs(value[i][i]) * sign;
            replace_row(j, value[i] * abs(value[j][i]) + value[j] * modified_value); // eliminates columns
        }
    }
    for (int i = 0; i < rows; ++i) // now we have a diagonal matrix with some numbers
                                   // this loop makes the values 1, making an identical matrix on the left
        if (value[i][i] != 0.0)
            replace_row(i, (value[i] / value[i][i]));

    // putting zero rows on the bottom
    for (int i = 0; i < rows - 1; ++i)
    {
        if (is_all_zero(i))
            swap_rows(i, i + 1);
    }
}
// returns A^(-1) from the form (I|A^(-1)
void Matrix::extract_inverse_matrix()
{

    int size = cols / 2;
    Matrix new_matrix(size, size);
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            new_matrix.value[i][j] = value[i][j + size];

    for (int i = 0; i < size; ++i)
        value[i].resize(size);

    set_cols(size);
    set_rows(size);

    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            value[i][j] = new_matrix.value[i][j];
}

void Matrix::extend_matrix_with_vector(vector<double> &b)
{
    set_cols(cols + 1);
    for (int i = 0; i < rows; ++i)
        value[i].resize(cols);

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
        {
            if (j == cols - 1)
            {
                value[i][j] = b[i];
            }
            else
            {
                value[i][j] = value[i][j];
            }
        }

    //
}

/*vector<vector<double>> Matrix::get_vec_inf_solutions() {
    vector<vector<double>> x(get_rows(), vector<double>(2, 0));

    for (int i = rows - 1; i > 0; --i) {
        if (bool(is_all_zero(i))) {
            x[i] = {0, 1};

        } else {
            // cout << "value[i][cols-1]: " << value[i][cols - 1] << '\t';
            vector<double> result{value[i][cols - 1], 0};

            for (int j = 0; j < cols - 1; ++j) {
                if (j == i) continue;
                result = result - value[i][j] * x[i];
            }
            // x[i] = result / value[i][i];
        }
        //
    }
    //

    return x;
}

// output
void Matrix::print_vec_inf_solutions() {
    vector<vector<double>> x(get_vec_inf_solutions());
    cout << "size" << x.size() << endl;
    for (int i = 0; i < x.size(); ++i) {
        if (x[i][1] == 0) {
            cout << x[i][0] << endl;
        } else {
            std ::string s = x[i][0] > 0 ? s = "+ " : s = "";

            cout << x[i][1] << "t " << s << x[i][0] << endl;
        }
    }
}
*/
bool Matrix::is_all_zero(int row_index)
{
    bool all_zero = true;
    for (int i = 0; i < cols; ++i)
    {
        if (value[row_index][i] != 0)
            return false;
    }

    return true;
}