#include<assemble_pressure_f.h>


void assemble_pressure_f_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                            Eigen::MatrixXd &M_particles, Eigen::VectorXd &f){
                                f = Eigen::VectorXd(M_particles.size());
                                //TODO: Implement B * q_j for the flattened grid
}