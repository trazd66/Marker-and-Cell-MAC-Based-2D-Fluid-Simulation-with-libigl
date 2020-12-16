#include<assemble_pressure_f.h>


void assemble_pressure_f_2d(double rho, double dx, double dy, double dt, 
                            Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                            Eigen::MatrixXd &M_signed_distance,Eigen::MatrixXd &M_particles,
                            Eigen::VectorXd &f){
                            assert(dx ==  dy && "x interval should equal to y interval.");
                            
                            f = Eigen::VectorXd((M_u.rows() - 1)*(M_v.cols() - 1));
                            Eigen::Matrix14d B;
                            B << -1, 1, -1, 1;      
                            Eigen::DiagonalMatrix<double,4> P_TP;
                            for (int i = 1; i < M_u.rows(); i++)
                            {
                                for (int j = 1; j < M_v.cols(); j++)
                                {
                                    if (M_signed_distance(i - 1,j) > 0)
                                        P_TP.diagonal()[0] = 0;
                                    if (M_signed_distance(i + 1,j) > 0)
                                        P_TP.diagonal()[1] = 0;
                                    if (M_signed_distance(i,j - 1) > 0)
                                        P_TP.diagonal()[2] = 0;
                                    if (M_signed_distance(i,j + 1) > 0)
                                        P_TP.diagonal()[3] = 0;
                                    Eigen::Vector4d q_j;
                                    q_j << M_u(i,j) , M_u(i + 1,j), M_v(i,j), M_v(i,j+1);
                                    f(i*M_particles.rows()+j) = B *P_TP* q_j;
                                }
                            }

                            f *= rho * dx / dt;
}