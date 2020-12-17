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

            f = Eigen::VectorXd(M_signed_distance.size());
            f.setZero();
            Eigen::Matrix14d B;
            B << -1, 1, -1, 1;      
            for (int x = 1; x < x_len_non_staggered; x++)
            {
                int x_non_staggered = x - 1;
                for (int y = 1; y < y_len_non_staggered; y++)
                {
                    double num_fluid_cells = 0.;
                    int y_non_staggered = y - 1;
                    int i_idx,j_idx;
                    
                    
                    get_matrix_index_2d(x_non_staggered,y_non_staggered,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
                    int idx = i_idx * x_len_non_staggered + j_idx;
                    Eigen::Vector4d q_j;
                    q_j.setZero();

                    get_matrix_index_2d(x_non_staggered,y-1,x_len_non_staggered,y_len,i_idx,j_idx);
                    if(M_signed_distance(i_idx,j_idx) < 0){//bottom
                        q_j[0] = M_v(i_idx,j_idx);
                    }

                    get_matrix_index_2d(x_non_staggered,y+1,x_len_non_staggered,y_len,i_idx,j_idx);
                    if(M_signed_distance(i_idx,j_idx) < 0){//top
                        q_j[1] = M_v(i_idx,j_idx);
                    }

                    get_matrix_index_2d(x+1,y_non_staggered,x_len,y_len_non_staggered,i_idx,j_idx);
                    if(M_signed_distance(i_idx,j_idx) < 0){//right
                        q_j[2] = M_u(i_idx,j_idx);
                    }

                    get_matrix_index_2d(x-1,y_non_staggered,x_len,y_len_non_staggered,i_idx,j_idx);
                    if(M_signed_distance(i_idx,j_idx) < 0){//left
                        q_j[3] = M_u(i_idx,j_idx);
                    }
                    
                    f[idx] = B * q_j;
                }
            }

            f *= rho * dx / dt;
}