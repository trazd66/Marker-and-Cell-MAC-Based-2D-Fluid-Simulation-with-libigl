#include <assemble_pressure_A.h>

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
            for (int i = 1; i < x_len_non_staggered; i++)
            {
                int i_non_staggered = i - 1;
                for (int j = 1; j < y_len_non_staggered; j++)
                {
                    double num_fluid_cells = 0.;
                    int j_non_staggered = j - 1;
                    int idx = i_non_staggered + j_non_staggered * y_len_non_staggered;//col-wise flattening

                    if(M_signed_distance(j-1,i_non_staggered) < 0){//top, saved at i th row and j -1 th col in Eigen
                        A_triplets.push_back(Eigen::Triplet<double>(idx,idx - y_len_non_staggered,-1.)); 
                        num_fluid_cells++;
                    }
                    if(M_signed_distance(j+1,i_non_staggered) < 0){//bottom
                        A_triplets.push_back(Eigen::Triplet<double>(idx,idx + y_len_non_staggered,-1.)); 
                        num_fluid_cells++;
                    }
                    if(M_signed_distance(j_non_staggered,i-1) < 0){//right
                        A_triplets.push_back(Eigen::Triplet<double>(idx,idx - 1,-1.)); 
                        num_fluid_cells++;
                    }
                    if(M_signed_distance(j_non_staggered,i+1) < 0){//left
                        A_triplets.push_back(Eigen::Triplet<double>(idx,idx + 1,-1.)); 
                        num_fluid_cells++;
                    }

                    A_triplets.push_back(Eigen::Triplet<double>(idx,idx,num_fluid_cells));
                }
            }
            A.setFromTriplets(A_triplets.begin(),A_triplets.end());
}