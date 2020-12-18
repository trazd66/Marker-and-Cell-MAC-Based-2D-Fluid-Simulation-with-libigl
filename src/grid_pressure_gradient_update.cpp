#include <grid_pressure_gradient_update.h>
#include <pressure_matrix_solve.h>
#include <utilities.h>
void grid_pressure_gradient_update_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v,
                                      Eigen::MatrixXd &M_particles,
                                      Eigen::MatrixXd &M_signed_distance,
                                      Eigen::SparseMatrixd &A,
                                      Eigen::VectorXd &f,
                                      double rho,
                                      double dt,
                                      double dx)
{

    Eigen::VectorXd p;
    solve_pressure_p(p, A, f);
    //TODO: calculate the pressure gradient using the pressure p
    //Construct big D to calculate the pressure gradient
    //then use the pressure gradient to update M_u and M_v

    int x_len = M_u.cols();
    int y_len = M_v.rows();
    int x_len_non_staggered = x_len - 1;
    int y_len_non_staggered = y_len - 1;

    int bb_size_x = M_v.rows();
    int bb_size_y = M_u.cols();

    int num_grids = x_len_non_staggered * y_len_non_staggered;
    Eigen::SparseMatrixd D(2 * num_grids, num_grids);
    std::vector<Eigen::Triplet<double>> D_triplets;

    for (int x = 0; x < x_len_non_staggered; x++)
    {
        for (int y = 0; y < y_len_non_staggered; y++)
        {
            int i_idx,j_idx;
            get_matrix_index_2d(x,y,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
            int idx = i_idx * y_len_non_staggered + j_idx;

            //v_x,y+1
            if(!on_boundary(y+1,y_len_non_staggered)){//top
                get_matrix_index_2d(x,y+1,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
                D_triplets.push_back(Eigen::Triplet<double>(2*idx,i_idx * y_len_non_staggered + j_idx,1.)); 
            }
            D_triplets.push_back(Eigen::Triplet<double>(2*idx,idx,-1.)); 

            //u_x+1,y
            if(!on_boundary(x+1,x_len_non_staggered)){//right
                get_matrix_index_2d(x+1,y,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
                D_triplets.push_back(Eigen::Triplet<double>(2*idx+1,i_idx * y_len_non_staggered + j_idx,1.)); 
            }
            D_triplets.push_back(Eigen::Triplet<double>(2*idx+1,idx,-1.)); 
            
        }
    }

    D.setFromTriplets(D_triplets.begin(), D_triplets.end());
    Eigen::VectorXd gradient_pressure(2 * num_grids);
    gradient_pressure = D * p;

    /* update M_u and M_v */

    std::cout << gradient_pressure << '\n';
    for (int x = 0; x < x_len_non_staggered; x++)
    {
        for (int y = 0; y < y_len_non_staggered; y++)
        {

            int i_idx,j_idx;
            get_matrix_index_2d(x,y,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
            int idx = i_idx * y_len_non_staggered + j_idx;

            //v_x,y+1
            get_matrix_index_2d(x,y+1,x_len_non_staggered,y_len,i_idx,j_idx);
            if(!on_boundary(y+1,y_len_non_staggered)){//top
                M_v(i_idx,j_idx) -= (dt / (rho * dx)) * gradient_pressure[2 * idx];
            }else{
                M_v(i_idx,j_idx) = 0;
            }

            //u_x+1,y
            get_matrix_index_2d(x+1,y,x_len,y_len_non_staggered,i_idx,j_idx);
            if(!on_boundary(x+1,x_len_non_staggered)){//right
                M_u(i_idx,j_idx) -= (dt / (rho * dx)) * gradient_pressure[2 * idx + 1];
            }else{
                M_u(i_idx,j_idx) = 0;

            }
        }
    }
    

}