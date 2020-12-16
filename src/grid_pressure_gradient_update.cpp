#include<grid_pressure_gradient_update.h>
#include<pressure_matrix_solve.h>


void grid_pressure_gradient_update_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v,
                                Eigen::MatrixXd &M_particles, 
                                Eigen::MatrixXd &M_signed_distance,
                                Eigen::SparseMatrixd &A,
                                Eigen::VectorXd &f){


                                Eigen::VectorXd p;
                                solve_pressure_p(p,A,f);
                                //TODO: calculate the pressure gradient using the pressure p
                                //Construct big D to calculate the pressure gradient
                                //then use the pressure gradient to update M_u and M_v

                                




                            



                                }