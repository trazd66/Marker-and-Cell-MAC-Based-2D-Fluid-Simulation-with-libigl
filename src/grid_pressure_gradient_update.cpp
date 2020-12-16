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

    int bb_size_x = M_v.rows();
    int bb_size_y = M_u.cols();
    int num_grids = bb_size_x * bb_size_y;
    Eigen::SparseMatrixd D(4 * num_grids, num_grids);
    std::vector<Eigen::Triplet<double>> D_triplets;

    for (int i = 0; i < bb_size_x; i++)
    {
        for (int j = 0; j < bb_size_y; j++)
        {
            /* p_x-1,y */
            if (i - 1 >= 0 && M_signed_distance(i - 1, j) < 0)
            {
                D_triplets.push_back(Eigen::Triplet<double>(4 * (bb_size_x * j + (i - 1)), bb_size_x * j + (i - 1), -1.0)); // p_x-1,y
                D_triplets.push_back(Eigen::Triplet<double>(4 * (bb_size_x * j + i), bb_size_x * j + i, 1.0));              // p_x,y
            }
            /* p_x+1,y */
            if (i + 1 < bb_size_x && M_signed_distance(i + 1, j) < 0)
            {
                D_triplets.push_back(Eigen::Triplet<double>(4 * (bb_size_x * j + (i + 1)), bb_size_x * j + (i + 1), 1.0)); // p_x+1,y
                D_triplets.push_back(Eigen::Triplet<double>(4 * (bb_size_x * j + i), bb_size_x * j + i, -1.0));            // p_x,y
            }
            /* p_x,y-1 */
            if (j - 1 >= 0 && M_signed_distance(i, j - 1) < 0)
            {
                D_triplets.push_back(Eigen::Triplet<double>(4 * (bb_size_x * j + i), bb_size_x * j + i, 1.0));              // p_x,y
                D_triplets.push_back(Eigen::Triplet<double>(4 * (bb_size_x * (j - 1) + i), bb_size_x * (j - 1) + i, -1.0)); // p_x,y-1
            }
            /* p_x,y+1 */
            if (j + 1 < bb_size_y && M_signed_distance(i, j + 1) < 0)
            {
                D_triplets.push_back(Eigen::Triplet<double>(4 * (bb_size_x * j + i), bb_size_x * j + i, -1.0));            // p_x,y
                D_triplets.push_back(Eigen::Triplet<double>(4 * (bb_size_x * (j + 1) + i), bb_size_x * (j + 1) + i, 1.0)); // p_x,y+1
            }
        }
    }

    D.setFromTriplets(D_triplets.begin(), D_triplets.end());
    Eigen::VectorXd gradient_pressure(4 * num_grids);
    gradient_pressure = D * p;

    /* update M_u and M_v */
    for (int i = 0; i < bb_size_x + 1; i++)
    {
        for (int j = 0; j < bb_size_y + 1; j++)
        {
            if (i == 0) {
                // left most grid, use p_x+1,y
                M_u(i, j) -= (dt / (rho * dx)) * gradient_pressure(4 * (bb_size_x * j + i) + 1);
            } else{
                // otherwise use gradient with respect to left grid
                M_u(i, j) -= (dt / (rho * dx)) * gradient_pressure(4 * (bb_size_x * j + i));
            }

            if (j == 0) {
                // bottom most grid, use p_x,y+1
                M_v(i, j) -= (dt / (rho * dx)) * gradient_pressure(4 * (bb_size_x * j + i) + 3);
            } else {
                // otherwise use gradient with respect to bottom grid
                M_v(i, j) -= (dt / (rho * dx)) * gradient_pressure(4 * (bb_size_x * j + i) + 2);
            }
        }
    }
}