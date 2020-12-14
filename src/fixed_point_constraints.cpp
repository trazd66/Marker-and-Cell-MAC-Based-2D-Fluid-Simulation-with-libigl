#include <fixed_point_constraints.h>
#include <algorithm>
#include <iostream>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    Eigen::SparseMatrixd PT= Eigen::SparseMatrixd();
    PT.resize(q_size,q_size - indices.size()*3);

    P = Eigen::SparseMatrixd(q_size - indices.size()*3,q_size);
    int j = 0;
    for (int i = 0; i < q_size/3; i++)
    {
	    if (std::count(indices.begin(), indices.end(), i) == 0){
            PT.coeffRef(3 * j,3 * j) = 1;
            PT.coeffRef(3 * j + 1,3 * j + 1) = 1;
            PT.coeffRef(3 * j + 2,3 * j + 2) = 1;
            j++;
        }
    }

    P = PT.transpose().eval();
}