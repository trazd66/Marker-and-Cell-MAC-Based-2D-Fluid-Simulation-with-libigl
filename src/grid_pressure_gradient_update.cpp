#include <grid_pressure_gradient_update.h>
#include <pressure_matrix_solve.h>

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

    int num_grids = M_signed_distance.size();
    Eigen::SparseMatrixd D(2 * num_grids, num_grids);
    std::vector<Eigen::Triplet<double>> D_triplets;

    for (int i = 1; i < x_len_non_staggered; i++)
    {
        int i_non_staggered = i - 1;
        for (int j = 1; j < y_len_non_staggered; j++)
        {
            int j_non_staggered = j - 1;
            int idx = i_non_staggered + j_non_staggered * y_len_non_staggered;//col-wise flattening

            if(M_signed_distance(j-1,i_non_staggered) < 0){//top, saved at i th row and j -1 th col in Eigen
                D_triplets.push_back(Eigen::Triplet<double>(4*idx,idx - y_len_non_staggered,1.)); 
                D_triplets.push_back(Eigen::Triplet<double>(4*idx,idx,-1.)); 
            }
            if(M_signed_distance(j_non_staggered,i-1) < 0){//right
                D_triplets.push_back(Eigen::Triplet<double>(4*idx+1,idx - 1,1.)); 
                D_triplets.push_back(Eigen::Triplet<double>(4*idx+1,idx,-1.)); 
            }
        }
    }

    D.setFromTriplets(D_triplets.begin(), D_triplets.end());
    Eigen::VectorXd gradient_pressure(2 * num_grids);
    gradient_pressure = D * p;

    /* update M_u and M_v */


    //distributing velocity to the staggerd grid M_u (w.r.t. x)
    for (int i = 1; i < x_len_non_staggered; i++)
    {
        int i_non_staggered = i - 1;
        for (int j = 1; j < y_len_non_staggered; j++)
        {
            int j_non_staggered = j - 1;
            int idx = i_non_staggered + j_non_staggered * y_len_non_staggered;//col-wise flattening

            M_u(j_non_staggered, j) -= (dt / (rho * dx)) * gradient_pressure(4 * idx);
            M_v(i, i_non_staggered) -= (dt / (rho * dx)) * gradient_pressure(4 * idx + 1);

        }
    }
}