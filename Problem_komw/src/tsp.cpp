#include "tsp.hpp"

double get_forbidden_cost() { return INF; }

double CostMatrix::min_in_row(int row, int forb_elem) {
    double min;

    for (int i = 0; i < int(costMatrix[row].size()); ++i) {
        if (i != forb_elem && !(std::isnan(costMatrix[row][i]))) {
            min = costMatrix[row][i];
            break;
        }
    }

    for (int i = 0; i < int(costMatrix[row].size()); ++i) {

        if (costMatrix[row][i] < min && i != forb_elem && !(std::isnan(costMatrix[row][i]))) {
            min = costMatrix[row][i];
        }
    }
    return min;
}

double CostMatrix::min_in_col(int col, int forb_elem) {
    double min;

    for (int i = 0; i < int(costMatrix.size()); ++i) {
        if (i != forb_elem && !(std::isnan(costMatrix[i][col]))) {
            min = costMatrix[i][col];
            break;
        }
    }

    for (int i = 0; i < int(costMatrix.size()); ++i) {
        if (costMatrix[i][col] < min && i != forb_elem && !(std::isnan(costMatrix[i][col]))) {
            min = costMatrix[i][col];
        }
    }
    return min;
}

void CostMatrix::reduce_rows() {
    std::vector<std::vector<double>> costMatrix_new;
    std::vector<double> new_row;
    for (int row = 0; row < int(costMatrix.size()); ++row) {
        double min = min_in_row(row);
        for (auto &e : costMatrix[row]) {
            e -= min;
            new_row.push_back(e);
        }
        costMatrix_new.push_back(new_row);
        new_row.clear();
    }
    costMatrix = costMatrix_new;
}

void CostMatrix::reduce_cols() {
    for (int col = 0; col < int(costMatrix.size()); ++col) {
        double min = min_in_col(col);
        for (int j = 0; j < int(size_col()); ++j) {
            costMatrix[j][col] -= min;
        }
    }
}

void CostMatrix::block_route(int row, int col) {
    costMatrix[row][col] = INF;
}

void CostMatrix::block_row(int row) {
    for (int j = 0; j < int(size_row()); ++j) {
        costMatrix[row][j] = INF;
    }
}

void CostMatrix::block_col(int col) {
    for (int j = 0; j < int(size_col()); ++j) {
        costMatrix[j][col] = INF;
    }
}

int CostMatrix::zero_in_row(int row, int forb_elem) {
    for (int j = 0; j < int(size_row()); ++j) {
        if (!(std::isnan(costMatrix[row][j])) && j != forb_elem) {
            return 1;
        }
    }
    return 0;
}

int CostMatrix::zero_in_col(int col, int forb_elem) {
    for (int j = 0; j < int(size_col()); ++j) {
        if (!(std::isnan(costMatrix[j][col])) && j != forb_elem) {
            return 1;
        }
    }
    return 0;
}

std::vector<int> tsp_step(CostMatrix costMatrix) {
    double max = 0;
    double e = 0;
    std::vector<int> vert;
    for (int i = 0; i < int(costMatrix.size_row()); ++i) {
        for (int j = 0; j < int(costMatrix.size_col()); ++j) {
            if (costMatrix[i][j] == 0) {
                e = costMatrix.min_in_row(i, j) + costMatrix.min_in_col(j, i);
                if (e > max) {
                    max = e;
                    vert.clear();
                    vert.push_back(i + 1);
                    vert.push_back(j + 1);
                }
            }
        }
    }
    return vert;
}

std::vector<std::vector<int>> tsp_2x2_matrix(CostMatrix costMatrix, std::vector<std::vector<int>> vert_dict) {
    std::vector<int> vert;

    for (int i = 0; i < int(costMatrix.size_row()); ++i) {
        for (int j = 0; j < int(costMatrix.size_col()); ++j) {
            if (costMatrix[i][j] == 0) {
                int a = costMatrix.zero_in_col(j, i);
                int b = costMatrix.zero_in_row(i, j);
                if (a + b == 1) {
                    vert.push_back(i + 1);
                    vert.push_back(j + 1);
                    vert_dict.push_back(vert);
                    vert.clear();
                }
            }
        }
    }
    return vert_dict;
}

std::vector<int> tsp_sort_vert(std::vector<std::vector<int>> vert_dict) {
    std::vector<int> vert_list(vert_dict[0]);

    for (int i = 0; i < int(vert_dict.size()) - 1; ++i) {
        for (auto e : vert_dict) {
            if (vert_list[vert_list.size() - 1] == e[0] && vert_list[vert_list.size() - 1] != vert_list[0]) {
                vert_list.push_back(e[1]);
            }
        }
    }
    return vert_list;
}

std::vector<int> tsp(std::vector<std::vector<double>> costMatrix_) {
    CostMatrix costMatrix(costMatrix_);

    std::vector<int> blocked_rows;
    std::vector<int> blocked_cols;

    std::vector<std::vector<int>> vert_dict;
    std::vector<int> vert;

    while (int((costMatrix.size_row()) - blocked_rows.size() > 2 && costMatrix.size_col() - blocked_cols.size() > 2)) {
        costMatrix.reduce_rows();
        costMatrix.reduce_cols();

        vert = tsp_step(costMatrix);

        costMatrix.block_route(vert[0] - 1, vert[1] - 1);
        costMatrix.block_route(vert[1] - 1, vert[0] - 1);

        costMatrix.block_row(vert[0] - 1);
        costMatrix.block_col(vert[1] - 1);

        vert_dict.push_back(vert);

        blocked_rows.push_back(vert[0]);
        blocked_cols.push_back(vert[1]);

        vert.clear();
    }

    costMatrix.reduce_rows();
    costMatrix.reduce_cols();

    sort(blocked_rows.begin(), blocked_rows.end());
    sort(blocked_cols.begin(), blocked_cols.end());

    vert_dict = tsp_2x2_matrix(costMatrix, vert_dict);

    std::vector<int> vert_list = tsp_sort_vert(vert_dict);

    return vert_list;
}



