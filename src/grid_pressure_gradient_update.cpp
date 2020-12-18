#include <grid_pressure_gradient_update.h>
#include <pressure_matrix_solve.h>
#include <utilities.h>
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
    //TODO: calculate the pressure gradient using the pressure p
    //Construct big D to calculate the pressure gradient
    //then use the pressure gradient to update M_u and M_v
    // std::cout << A.norm() << '\n';
    // std::cout << f.norm() << '\n';
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
            int i_idx,j_idx;
            get_matrix_index_2d(x,y,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
            int idx = i_idx * y_len_non_staggered + j_idx;
            M_pressure(i_idx,j_idx) = p[idx];
        }
    }

    for (int x = 0; x < x_len_non_staggered; x++)
    {
        for (int y = 0; y < y_len_non_staggered; y++)
        {
            int i_idx,j_idx;
            get_matrix_index_2d(x,y,x_len_non_staggered,y_len_non_staggered,i_idx,j_idx);
			if (i_idx != 0) {

				double du = (dt / (rho * dx)) * (M_pressure(i_idx,j_idx) - M_pressure(i_idx-1,j_idx));
				M_u(i_idx,j_idx)-= du;
			}
			if (j_idx != 0) {

				double dv = (dt / (rho * dx)) * (M_pressure(i_idx,j_idx) - M_pressure(i_idx,j_idx-1));
				M_v(i_idx,j_idx) -= dv;
			}
        }
    }

   

}