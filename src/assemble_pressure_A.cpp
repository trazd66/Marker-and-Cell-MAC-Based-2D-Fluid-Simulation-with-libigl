#include <assemble_pressure_A.h>

void assemble_pressure_A_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                            Eigen::MatrixXd &M_particles, 
                            Eigen::MatrixXd &M_signed_distance,
                            Eigen::SparseMatrixd &A){

            A = Eigen::SparseMatrixd(M_u.rows() - 1,M_v.cols() - 1);                            
            std::vector<Eigen::Triplet<double>> A_triplets;

            int row_len = M_particles.rows();
            int col_len = M_particles.cols();
            for (int i = 1; i < M_u.rows(); i++)
            {
                for (int j = 1; j < M_v.rows(); j++)
                {
                    double num_fluid_cells = 0.;
                    int idx = (i-1)*row_len + (j-1);
                    if(M_signed_distance(i,j+1) < 0){//top
                        A_triplets.push_back(Eigen::Triplet<double>(idx,i*row_len+(j+1),-1.)); 
                        num_fluid_cells++;
                    }
                    if(M_signed_distance(i,j-1) < 0){//bottom
                        A_triplets.push_back(Eigen::Triplet<double>(idx,i*row_len+(j-1),-1.)); 
                        num_fluid_cells++;
                    }
                    if(M_signed_distance(i+1,j) < 0){//right
                        A_triplets.push_back(Eigen::Triplet<double>(idx,(i+1)*row_len+j,-1.)); 
                        num_fluid_cells++;
                    }
                    if(M_signed_distance(i-1,j) < 0){//left
                        A_triplets.push_back(Eigen::Triplet<double>(idx,(i-1)*row_len+j,-1.)); 
                        num_fluid_cells++;
                    }

                    A_triplets.push_back(Eigen::Triplet<double>(idx,idx,num_fluid_cells));
                }
            }
            A.setFromTriplets(A_triplets.begin(),A_triplets.end());
}