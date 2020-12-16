#ifndef EIGENTYPES_H
#define EIGENTYPES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace Eigen {

    //dense types
    using Vector6d = Eigen::Matrix<double, 6,1>;
    using Matrix15d = Eigen::Matrix<double,1,5>;
    using Matrix14d = Eigen::Matrix<double,1,4>;
    using Matrix36d = Eigen::Matrix<double, 3,6>;
    using Matrix66d  = Eigen::Matrix<double, 6,6>;
    using Matrix44f = Eigen::Matrix<float, 4,4>;
    using Matrix45d = Eigen::Matrix<double,4,5>;
    //sparse types
    using SparseMatrixd = Eigen::SparseMatrix<double>;

    
    //Tensor types
    struct TensorXd{
        std::vector<Eigen::MatrixXd> data;

        void init_tensor(const int x,const int y,const int z){
            for (size_t i = 0; i < x; i++)
            {
                data.push_back(Eigen::MatrixXd(y,z));
                data[i].setZero();
            }            
        };
    };

}

#endif 
