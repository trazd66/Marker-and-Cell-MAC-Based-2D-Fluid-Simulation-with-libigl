#include <grid_pressure_gradient_update.h>
#include <pressure_matrix_solve.h>
#include <utilities.h>

/***
 * Calculates grid pressures using the equation from lecture slides.
 * Then updates grid velocities using the new pressures.
 * https://github.com/dilevin/CSC417-physics-based-animation/blob/master/lectures/10-fluid-simulation-final.pdf
 * Effect: updates M_u and M_v.
 ***/
void grid_pressure_gradient_update_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v,
                                      Eigen::MatrixXd &M_particles,
                                      Eigen::MatrixXd &M_pressure,
                                      Eigen::MatrixXd &M_fluid,
                                      Eigen::SparseMatrixd &A,
                                      Eigen::VectorXd &f,
                                      double rho,
                                      double dt,
                                      double dx)
{

    Eigen::VectorXd p;
    solve_pressure_p(p, A, f);

    // Calculate the pressure gradient using the pressure p
    // Construct big D to calculate the pressure gradient
    // then use the pressure gradient to update M_u and M_v
    int x_len = M_u.cols();
    int y_len = M_v.rows();
    int x_len_non_staggered = x_len - 1;
    int y_len_non_staggered = y_len - 1;

    int bb_size_x = M_v.rows();
    int bb_size_y = M_u.cols();


    for (int x = 0; x < x_len_non_staggered; x++)
    {
        for (int y = 0; y < y_len_non_staggered; y++)
        {
            int i_idx, j_idx;
            get_matrix_index_2d(x, y, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
            int idx = i_idx * y_len_non_staggered + j_idx;
            M_pressure(i_idx, j_idx) = p[idx];
        }
    }

    // updates grid velocity using the new pressure
    for (int x = 0; x < x_len_non_staggered; x++)
    {
        for (int y = 0; y < y_len_non_staggered; y++)
        {
            int i_idx, j_idx, i,j;
            get_matrix_index_2d(x, y, x_len_non_staggered, y_len_non_staggered, i,j);
            if (!on_boundary(x-1,x_len_non_staggered))
            {
                get_matrix_index_2d(x - 1, y, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                double du = (dt / (rho * dx)) * ( M_pressure(i,j) - M_pressure(i_idx, j_idx));
                M_u(i,j) -= du;
            }
            if (!on_boundary(y-1,y_len_non_staggered))
            {
                get_matrix_index_2d(x , y - 1, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                double dv = (dt / (rho * dx)) * (M_pressure(i,j) - M_pressure(i_idx, j_idx));
                M_v(i,j) -= dv;
            }
        }
    }
}
