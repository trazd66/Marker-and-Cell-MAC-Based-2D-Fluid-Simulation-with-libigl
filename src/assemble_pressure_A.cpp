#include <assemble_pressure_A.h>

void assemble_pressure_A_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                            Eigen::MatrixXd &M_particles, Eigen::SparseMatrixd &A){

                                A = Eigen::SparseMatrixd(M_particles.size(),5);
                                //TODO : Implement the per element B_j ,D_j matrix
}