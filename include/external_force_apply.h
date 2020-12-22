#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    Apply external forces using the time step value, currently its only gravity.
    Effect: updates M_u and M_v.
*/
void apply_external_force_2d(Eigen::MatrixXd &M_u,Eigen::MatrixXd &M_v, double dt);
