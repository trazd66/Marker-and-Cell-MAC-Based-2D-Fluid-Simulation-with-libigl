#include<assemble_pressure_f.h>
#include<utilities.h>

void assemble_pressure_f_2d(double rho, double dx, double dy, double dt, 
                            Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                            Eigen::MatrixXd &M_signed_distance,Eigen::MatrixXd &M_particles,
                            Eigen::VectorXd &f){
            assert(dx ==  dy && "x interval should equal to y interval.");
            int x_len = M_u.cols();
            int y_len = M_v.rows();
            int x_len_non_staggered = x_len - 1;
            int y_len_non_staggered = y_len - 1;

            f = Eigen::VectorXd(x_len_non_staggered * y_len_non_staggered);
            f.setZero();
            Eigen::Matrix14d B;
            B << -1, 1, -1, 1;      
            for (int x = 0; x < x_len_non_staggered; x++)
            {
                for (int y = 0; y < y_len_non_staggered; y++)
                {

                    int i_idx,j_idx;
                    get_matrix_index_2d(x,y,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
                    int idx = i_idx * y_len_non_staggered + j_idx;

                    Eigen::Vector4d q_j;
                    q_j.setZero();

                    //v_xy
                    if(!on_boundary(y,y_len)){//bottom
                        get_matrix_index_2d(x,y,x_len_non_staggered,y_len,i_idx,j_idx);
                        q_j[0] = M_v(i_idx,j_idx);
                    }

                    //v_x,y+1
                    if(!on_boundary(y+1,y_len)){//top
                        get_matrix_index_2d(x,y+1,x_len_non_staggered,y_len,i_idx,j_idx);
                        q_j[1] = M_v(i_idx,j_idx);
                    }

                    //u_x,y
                    if(!on_boundary(x,x_len)){//left
                        get_matrix_index_2d(x,y,x_len,y_len_non_staggered,i_idx,j_idx);
                        q_j[2] = M_u(i_idx,j_idx);
                    }

                    //u_x+1,y
                    if(!on_boundary(x+1,x_len)){//right
                        get_matrix_index_2d(x+1,y,x_len,y_len_non_staggered,i_idx,j_idx);
                        q_j[3] = M_u(i_idx,j_idx);
                    }

                    f[idx] = B * q_j;
                }
            }

            f *= -rho * dx / dt;

}