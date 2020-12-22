#include <update_markers.h>
#include <utilities.h>
#include <iostream>

/***
 *
 * Given the u and v velocity grid, update the M_liquid matrix from all the marker particles
 * Effect: updates M_fluid
 ***/
void update_markers_2d(Eigen::MatrixXd &M_particles, double grid_interval,
                       Eigen::MatrixXd &M_fluid)
{

    M_fluid.setZero();
    for (int i = 0; i < M_particles.rows(); i++)
    {
        int i_idx, j_idx;
        get_matrix_index_2d((int)(M_particles(i, 0) / grid_interval), (int)(M_particles(i, 1) / grid_interval), M_fluid.cols(), M_fluid.rows(), i_idx, j_idx);
        ASSERT((i_idx >= 0) && (i_idx <= M_fluid.rows() - 1), "i_idx -> " << i_idx);
        ASSERT((j_idx >= 0) && (j_idx <= M_fluid.cols() - 1), "j_idx -> " << j_idx);

        M_fluid(i_idx, j_idx) = 1;
    }
}
