#ifndef SOURCE_FILES
#define SOURCE_FILES

#include <cmath>
#include <functional>
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <algorithm>

#define INF (NAN)


std::vector<int> tsp(std::vector<std::vector<double>> costMatrix_);

double get_forbidden_cost();

class CostMatrix {
private:
    std::vector<std::vector<double>> costMatrix;

public:
    CostMatrix(std::vector<std::vector<double>> const &cost_matrix):  costMatrix(cost_matrix) {};

    const std::vector<double> &operator[](std::size_t index) const { return costMatrix[index]; }

    std::vector<double> &operator[](std::size_t index) { return costMatrix[index]; }

    double min_in_row(int row, int forb_elem = -1);

    double min_in_col(int col, int forb_elem = -1);

    void reduce_rows();

    void reduce_cols();

    void block_route(int row, int col);

    void block_row(int row);

    void block_col(int col);

    int zero_in_row(int row, int forb_elem = -1);

    int zero_in_col(int col, int forb_elem = -1);

    std::size_t size_row() const { return costMatrix.size(); }

    std::size_t size_col() const { return costMatrix[0].size(); }
};

std::vector<int> tsp_step(CostMatrix costMatrix);

std::vector<std::vector<int>> tsp_2x2_matrix(CostMatrix costMatrix, std::vector<std::vector<int>> vert_dict);

std::vector<int> tsp_sort_vert(std::vector<std::vector<int>> vert_dict);

#endif