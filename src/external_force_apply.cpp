#include <external_force_apply.h>

// gravity acceleration
Eigen::Vector2d g(0., -9.8);

/*
    Apply external forces using the time step value, currently its only gravity.
    Effect: updates M_u and M_v.
*/
void apply_external_force_2d(Eigen::MatrixXd &M_u,Eigen::MatrixXd &M_v, double dt){
    //gravity
    M_v += Eigen::MatrixXd::Ones(M_v.rows(),M_v.cols()) * g[1] * dt; 
}
