#include <Eigen/Sparse>
#include <EigenTypes.h>

/***
 * 
 * 
 * 
 * 
 ***/
void grid_pressure_gradient_update_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                                Eigen::MatrixXd &M_particles, 
                                Eigen::MatrixXd &M_signed_distance,
                                Eigen::SparseMatrixd &A,
                                Eigen::VectorXd &f);