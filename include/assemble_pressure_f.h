#include <Eigen/Sparse>
#include <EigenTypes.h>

/***
 * Assembles the global pressure vector f in the 2d case, where f_j = (rho / dt) * B *(P^TP) * q_j
 *
 * input:
 * double rho
 * double dx
 * double dy
 * double dt
 * Matrix (bb_size_x X bb_size_y + 1 ) M_u
 * Matrix (bb_size_x + 1 X bb_size_y ) M_v
 * Matrix (bb_size_x X bb_size_y) M_fluid marker matrix that marks fluid cells and air cells
 * MatrixXd (num_particles X 3) M_particles
 *
 * output:
 * VectorXd (bb_size_x * bb_size_y) f  //flattened grid
 ***/
void assemble_pressure_f_2d(double rho, double dx, double dy, double dt, 
                            Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                            Eigen::MatrixXd &M_fluid,Eigen::MatrixXd &M_particles,
                            Eigen::VectorXd &f);
