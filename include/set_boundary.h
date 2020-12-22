#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    Update grid velocity such that boundary grids have velocity of zero.
    Effect: updates M_u and M_v.
*/
void set_boundary_2d(Eigen::MatrixXd &M_u,Eigen::MatrixXd &M_v);
