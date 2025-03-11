#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED
#include <vector>
using std::vector;

class Matrix {
public:
    vector<vector<double>> value;
    Matrix(int rows, int cols) : rows(rows), cols(cols), value(rows, vector<double>(cols)) {}

    Matrix(vector<vector<double>>& other) : value(other) {}

    Matrix(Matrix& other) {
        this->value = other.value;
    }

    void matrix_from_vector(vector<vector<double>>& vec) {
        value = vec;
    }

    int get_rows();
    void set_rows(int n);

    int get_cols();
    void set_cols(int n);

    void print();

    vector<vector<double>> get_minor(int deleted_row, int deleted_column);

    double det();
    void swap_rows(int a, int b);

    vector<double> add_rows_by_index(int a, int b);

    void add_identical_matrix();
    void replace_row(int row_number, vector<double> replacement_row);

    void gauss_jordan();
    void extract_inverse_matrix();

    void extend_matrix_with_vector(vector<double>& b);

    vector<vector<double>> to_vector() { return value; }

    vector<vector<double>> get_vec_inf_solutions();

    void print_vec_inf_solutions();
    bool is_all_zero(int row_index);

private:
    int rows;
    int cols;
};
#endif
