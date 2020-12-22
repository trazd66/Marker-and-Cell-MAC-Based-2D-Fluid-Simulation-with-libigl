#include <external_force_apply.h>


Eigen::Vector2d g(0., -9.8); // gravity acceleration

void apply_external_force_2d(Eigen::MatrixXd &M_u,Eigen::MatrixXd &M_v, double dt){
    //gravity
    M_v += Eigen::MatrixXd::Ones(M_v.rows(),M_v.cols()) * g[1] * dt; 
}