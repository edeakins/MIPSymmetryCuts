#ifndef LINEAR_PROGRAM_H_
#define LINEAR_PROGRAM_H_

#include <vector>

class linear_program{
public:
    int num_col;
    int num_row;
    int num_total;
    std::vector<double> a_value;
    std::vector<int> a_index;
    std::vector<int> a_start;
    std::vector<double> ar_value;
    std::vector<int> ar_index;
    std::vector<int> ar_start;
    std::vector<double> col_lower;
    std::vector<double> col_upper;
    std::vector<double> row_lower;
    std::vector<double> row_upper;
    std::vector<double> col_cost;
};

#endif