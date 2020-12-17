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

            A = Eigen::SparseMatrixd(M_signed_distance.size(),M_signed_distance.size());//non-staggered                            
            std::vector<Eigen::Triplet<double>> A_triplets;

            //note that although in Eigen::Matrix form, M(i,j) is i th row, j th column
            //but in the paper, M(i,j) is i th column, j th row, which we will use here 
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
                    // std::cout <<"i_idx: " << i_idx <<" j_idx: " << j_idx << '\n';

                    get_matrix_index_2d(x_non_staggered,y+1,x_len_non_staggered,y_len,i_idx,j_idx);
                    if(M_signed_distance(i_idx,j_idx) < 0){//top
                        A_triplets.push_back(Eigen::Triplet<double>(idx,i_idx * x_len_non_staggered + j_idx,-1.)); 
                        num_fluid_cells++;
                    }

                    get_matrix_index_2d(x_non_staggered,y-1,x_len_non_staggered,y_len,i_idx,j_idx);
                    if(M_signed_distance(i_idx,j_idx) < 0){//bottom
                        A_triplets.push_back(Eigen::Triplet<double>(idx,i_idx * x_len_non_staggered + j_idx,-1.)); 
                        num_fluid_cells++;
                    }

                    get_matrix_index_2d(x+1,y_non_staggered,x_len,y_len_non_staggered,i_idx,j_idx);
                    if(M_signed_distance(i_idx,j_idx) < 0){//right
                        A_triplets.push_back(Eigen::Triplet<double>(idx,i_idx * x_len_non_staggered + j_idx,-1.)); 
                        num_fluid_cells++;
                    }

                    get_matrix_index_2d(x-1,y_non_staggered,x_len,y_len_non_staggered,i_idx,j_idx);
                    if(M_signed_distance(i_idx,j_idx) < 0){//left
                        A_triplets.push_back(Eigen::Triplet<double>(idx,i_idx * x_len_non_staggered + j_idx,-1.)); 
                        num_fluid_cells++;
                    }

                    A_triplets.push_back(Eigen::Triplet<double>(idx,idx,num_fluid_cells));
                }
            }
            A.setFromTriplets(A_triplets.begin(),A_triplets.end());
}