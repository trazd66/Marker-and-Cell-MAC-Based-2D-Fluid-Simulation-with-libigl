#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <iostream>

/*
    calculate grid velocity extrapolation, transfers fluid cell velocity onto air cells.
    Effect: updates M_u and M_v.
*/
void extrapolate_velocity_2d(Eigen::MatrixXd & M_u,Eigen::MatrixXd & M_v, Eigen::MatrixXd & M_fluid);
