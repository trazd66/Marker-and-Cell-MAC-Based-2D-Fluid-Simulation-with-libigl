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
                                      Eigen::MatrixXd &M_pressure,
                                      Eigen::MatrixXd &M_fluid,
                                      Eigen::SparseMatrixd &A,
                                      Eigen::VectorXd &f,
                                      double rho,
                                      double dt,
                                      double dx);