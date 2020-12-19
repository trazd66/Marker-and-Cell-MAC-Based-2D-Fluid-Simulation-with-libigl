#include <Eigen/Sparse>
#include <EigenTypes.h>


void advect_particle_2d(Eigen::MatrixXd &M_particles, Eigen::VectorXd &u, Eigen::VectorXd &v, double dt, Eigen::MatrixXd M_u, Eigen::MatrixXd M_v, double grid_interval, int len_x, int len_y);