#include<assemble_pressure_f.h>


void assemble_pressure_f_2d(double rho, double dx, double dy, double dt, 
                            Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                            Eigen::MatrixXd &M_signed_distance,Eigen::MatrixXd &M_particles,
                            Eigen::VectorXd &f){
                            assert(dx ==  dy && "x interval should equal to y interval.");
                            int x_len = M_u.cols();
                            int y_len = M_v.rows();
                            int x_len_non_staggered = x_len - 1;
                            int y_len_non_staggered = y_len - 1;

                            f = Eigen::VectorXd(M_signed_distance.size());
                            Eigen::Matrix14d B;
                            B << -1, 1, -1, 1;      
                            Eigen::DiagonalMatrix<double,4> P_TP;
                            for (int i = 1; i < x_len_non_staggered; i++)
                            {
                                int i_non_staggered = i - 1;
                                for (int j = 1; j < y_len_non_staggered; j++)
                                {
                                    int j_non_staggered = j - 1;
                                    int idx = i_non_staggered + j_non_staggered * y_len_non_staggered;//col-wise flattening

                                    if (M_signed_distance(j_non_staggered,i - 1) > 0)
                                        P_TP.diagonal()[0] = 0;
                                    if (M_signed_distance(j_non_staggered,i + 1) > 0)
                                        P_TP.diagonal()[1] = 0;
                                    if (M_signed_distance(j - 1,i_non_staggered) > 0)
                                        P_TP.diagonal()[2] = 0;
                                    if (M_signed_distance(j + 1,i_non_staggered) > 0)
                                        P_TP.diagonal()[3] = 0;
                                    Eigen::Vector4d q_j;
                                    q_j << M_u(j_non_staggered,i) , M_u(j_non_staggered,i + 1), M_v(j,i_non_staggered), M_v(j+1,i_non_staggered);
                                    f(idx) = B *P_TP* q_j;
                                }
                            }

                            f *= rho * dx / dt;
}