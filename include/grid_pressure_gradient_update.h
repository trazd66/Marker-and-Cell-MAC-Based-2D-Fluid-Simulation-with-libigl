#include <Eigen/Sparse>
#include <EigenTypes.h>

/***
 * Calculates grid pressures using the equation from lecture slides.
 * Then updates grid velocities using the new pressures.
 * https://github.com/dilevin/CSC417-physics-based-animation/blob/master/lectures/10-fluid-simulation-final.pdf
 * Effect: updates M_u and M_v.
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
