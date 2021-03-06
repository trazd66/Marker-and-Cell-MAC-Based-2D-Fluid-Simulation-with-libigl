#include <assemble_pressure_A.h>
#include <utilities.h>

/***
 * Assembles the global pressure matrix A in the 2d case, where A_j = B*(P^TP)*D
 *
 * input:
 * Matrix (bb_size_x X bb_size_y + 1 ) M_u
 * Matrix (bb_size_x + 1 X bb_size_y ) M_v
 * MatrixXd (num_particles X 3) M_particles
 * Matrix (bb_size_x X bb_size_y) M_fluid marker matrix that marks fluid cells and air cells
 *
 * output:
 * SparseMatrixd (bb_size_x * bb_size_y X 5) A  //flattened grid
 ***/
void assemble_pressure_A_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v,
                            Eigen::MatrixXd &M_particles,
                            Eigen::MatrixXd &M_fluid,
                            Eigen::SparseMatrixd &A)
{

    int x_len = M_u.cols();
    int y_len = M_v.rows();
    int x_len_non_staggered = x_len - 1;
    int y_len_non_staggered = y_len - 1;

    A = Eigen::SparseMatrixd(y_len_non_staggered * x_len_non_staggered, y_len_non_staggered * x_len_non_staggered); //non-staggered
    std::vector<Eigen::Triplet<double>> A_triplets;

    //note that although in Eigen::Matrix form, M(i,j) is i th row, j th column
    //but in the paper, M(i,j) is i th column, j th row, which we will use here
    for (int x = 0; x < x_len_non_staggered; x++)
    {
        for (int y = 0; y < y_len_non_staggered; y++)
        {
            double num_fluid_cells = 0.;

            int i_idx, j_idx, center_idx_i, center_idx_j;
            get_matrix_index_2d(x, y, x_len_non_staggered, y_len_non_staggered, center_idx_i, center_idx_j);
            int idx = center_idx_i * y_len_non_staggered + center_idx_j;
            //v_x,y+1
            if (!on_boundary(y + 1, y_len_non_staggered))
            { //top
                num_fluid_cells++;
            }
            //v_x,y-1
            if (!on_boundary(y - 1, y_len_non_staggered))
            { //bot
                num_fluid_cells++;
            }
            //u_x-1,y
            if (!on_boundary(x - 1, x_len_non_staggered))
            { //left
                num_fluid_cells++;
            }
            //u_x+1,y
            if (!on_boundary(x + 1, x_len_non_staggered))
            { //right
                num_fluid_cells++;
            }
            A_triplets.push_back(Eigen::Triplet<double>(idx, idx, num_fluid_cells));

            if (M_fluid(center_idx_i, center_idx_j) == 1)
            {

                get_matrix_index_2d(x + 1, y, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (!on_boundary(x + 1, x_len_non_staggered) && M_fluid(i_idx, j_idx) == 1)
                {

                    A_triplets.push_back(Eigen::Triplet<double>(idx, i_idx * y_len_non_staggered + j_idx, -1.));
                    A_triplets.push_back(Eigen::Triplet<double>(i_idx * y_len_non_staggered + j_idx, idx, -1.)); //symmetric
                }

                get_matrix_index_2d(x, y + 1, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (!on_boundary(y + 1, y_len_non_staggered) && M_fluid(i_idx, j_idx) == 1)
                {
                    A_triplets.push_back(Eigen::Triplet<double>(idx, i_idx * y_len_non_staggered + j_idx, -1.));
                    A_triplets.push_back(Eigen::Triplet<double>(i_idx * y_len_non_staggered + j_idx, idx, -1.)); //symmetric
                }
            }
        }
    }
    A.setFromTriplets(A_triplets.begin(), A_triplets.end());

}
