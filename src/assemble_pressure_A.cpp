#include <assemble_pressure_A.h>
#include <utilities.h>

void assemble_pressure_A_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                            Eigen::MatrixXd &M_particles, 
                            Eigen::MatrixXd &M_signed_distance,
                            Eigen::SparseMatrixd &A){

        int x_len = M_u.cols();
        int y_len = M_v.rows();
        int x_len_non_staggered = x_len - 1;
        int y_len_non_staggered = y_len - 1;

        A = Eigen::SparseMatrixd(y_len_non_staggered * x_len_non_staggered,y_len_non_staggered * x_len_non_staggered);//non-staggered                            
        std::vector<Eigen::Triplet<double>> A_triplets;

        //note that although in Eigen::Matrix form, M(i,j) is i th row, j th column
        //but in the paper, M(i,j) is i th column, j th row, which we will use here 
        for (int x = 0; x < x_len_non_staggered; x++)
        {
            for (int y = 0; y < y_len_non_staggered; y++)
            {
                double num_fluid_cells = 0.;

                int i_idx,j_idx;
                get_matrix_index_2d(x,y,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
                int idx = i_idx * y_len_non_staggered + j_idx;

                //v_x,y+1
                if(!on_boundary(y+1,y_len_non_staggered)){//top
                    get_matrix_index_2d(x,y+1,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
                    A_triplets.push_back(Eigen::Triplet<double>(idx,i_idx * y_len_non_staggered + j_idx,-1.)); 
                    num_fluid_cells++;
                }
                //v_x,y-1
                if(!on_boundary(y-1,y_len_non_staggered)){//bot
                    get_matrix_index_2d(x,y-1,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
                    A_triplets.push_back(Eigen::Triplet<double>(idx,i_idx * y_len_non_staggered + j_idx,-1.)); 
                    num_fluid_cells++;
                }

                //u_x-1,y
                if(!on_boundary(x-1,x_len_non_staggered)){//left
                    get_matrix_index_2d(x-1,y,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
                    A_triplets.push_back(Eigen::Triplet<double>(idx,i_idx * y_len_non_staggered + j_idx,-1.)); 
                    num_fluid_cells++;
                }
                //u_x+1,y
                if(!on_boundary(x+1,x_len_non_staggered)){//right
                    get_matrix_index_2d(x+1,y,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
                    A_triplets.push_back(Eigen::Triplet<double>(idx,i_idx * y_len_non_staggered + j_idx,-1.)); 
                    num_fluid_cells++;
                }
                
                A_triplets.push_back(Eigen::Triplet<double>(idx,idx,num_fluid_cells));
            }
        }

        // std::cout << A << '\n';
        A.setFromTriplets(A_triplets.begin(),A_triplets.end());
}