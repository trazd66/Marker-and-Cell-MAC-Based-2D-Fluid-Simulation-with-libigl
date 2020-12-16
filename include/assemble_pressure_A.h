#include <Eigen/Sparse>
#include <EigenTypes.h>

/***
 * Assembles the global pressure matrix A in the 2d case, where A_j = B*(P^TP)*D
 * 
 * input:
 * Matrix (bb_size_x X bb_size_y + 1 ) M_u
 * Matrix (bb_size_x + 1 X bb_size_y ) M_v
 * MatrixXd (num_particles X 3) M_particles
 * Matrix (bb_size_x X bb_size_y) M_signed_distance, the signed distance matrix for the grid
 *
 * output:
 * SparseMatrixd (bb_size_x * bb_size_y X 5) A  //flattened grid
 ***/
void assemble_pressure_A_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                            Eigen::MatrixXd &M_particles, 
                            Eigen::MatrixXd &M_signed_distance,
                            Eigen::SparseMatrixd &A);