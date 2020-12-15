#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/***
 * input:
 * double bb_size_x  x dimension of the bounding box
 * double bb_size_y  y dimension of the bounding box
 * double bb_size_z  z dimension of the bounding box
 * double grid_interval size of the grid interval, determines the number of grid cells
 * double num_particles number of particles
 * 
 * output:
 * Matrix (num_particles X 3) M_particles the particle matrix 
 * Tensor (bb_size_x X bb_size_y + 1 X bb_size_z) M_u a 3D Tensor that contains the x velociies of the grid
 * Tensor (bb_size_x X bb_size_y X bb_size_z + 1) M_v a 3D Tensor that contains the y velociies of the grid
 * Tensor (bb_size_x + 1 X bb_size_y X bb_size_z) M_w a 3D Tensor that contains the z velociies of the grid
 ***/
void init_state_3d(const int bb_size_x,const int bb_size_y,const int bb_size_z, 
                const double grid_interval, const int num_particles, Eigen::MatrixXd &M_particles,
                Eigen::TensorXd &M_u, Eigen::TensorXd &M_v, Eigen::TensorXd &M_w);

/***
 * input:
 * double bb_size_x  x dimension of the bounding box
 * double bb_size_y  y dimension of the bounding box
 * double grid_interval size of the grid interval, determines the number of grid cells
 * double num_particles number of particles
 * 
 * output:
 * Matrix (num_particles X 2) the particle matrix 
 * Matrix (bb_size_x X bb_size_y + 1 ) M_u a 2D matrix that contains the x velociies of the grid
 * Matrix (bb_size_x + 1 X bb_size_y ) M_v a 2D matrix that contains the y velociies of the grid
 ***/
void init_state_2d(const int bb_size_x,const int bb_size_y, 
                const double grid_interval, const int num_particles, Eigen::MatrixXd &M_particles,
                Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v, Eigen::MatrixXd &M_pressures);