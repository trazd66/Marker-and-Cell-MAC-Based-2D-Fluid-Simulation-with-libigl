#include <particle_to_grid.h>
#include <iostream>
#include <utilities.h>
/*
    Accumulate u of the single particle onto the whole grid matrix using bilinear interpolation.
    There are four nodes affected by the particle, four corners correspondingly.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    EFFECT: Updates M_u
*/
void u_particle_onto_grid_u(Eigen::MatrixXd &M_u, Eigen::Vector2d &pos_particle, double u_particle, double grid_interval, const int len_x,const int len_y) {

    if (is_out_of_boundary(pos_particle, grid_interval, len_x, len_y)) return;

    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)(particle_x / grid_interval);
    int grid_y_start = (int)((particle_y - 0.5 * grid_interval) / grid_interval);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    Eigen::Vector4d bilinear_coeff;
    get_bilinear_coeff(bilinear_coeff, grid_interval, grid_x_start, grid_x_end, grid_y_start, grid_y_end, particle_x, particle_y);

    int i_idx, j_idx;
    get_matrix_index_2d(grid_x_start, grid_y_start, len_x+1, len_y, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y - 1) && (j_idx >= 0) && (j_idx <= len_x))
        M_u(i_idx, j_idx) += bilinear_coeff[0] * u_particle;

    get_matrix_index_2d(grid_x_start, grid_y_end, len_x+1, len_y, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y - 1) && (j_idx >= 0) && (j_idx <= len_x))
        M_u(i_idx, j_idx) += bilinear_coeff[1] * u_particle;

    get_matrix_index_2d(grid_x_end, grid_y_start, len_x+1, len_y, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y - 1) && (j_idx >= 0) && (j_idx <= len_x))
        M_u(i_idx, j_idx) += bilinear_coeff[2] * u_particle;

    get_matrix_index_2d(grid_x_end, grid_y_end, len_x+1, len_y, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y - 1) && (j_idx >= 0) && (j_idx <= len_x))
        M_u(i_idx, j_idx) += bilinear_coeff[3] * u_particle;
}

/*
    Accumulate v of the single particle onto the whole grid matrix using bilinear interpolation.
    There are four nodes affected by the particle, four corners correspondingly.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    EFFECT: Updates M_v
*/
void v_particle_onto_grid_v(Eigen::MatrixXd &M_v, Eigen::Vector2d &pos_particle, double v_particle, double grid_interval, const int len_x,const int len_y)
{

    if (is_out_of_boundary(pos_particle, grid_interval, len_x, len_y)) return;

    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)((particle_x - 0.5 * grid_interval) / grid_interval);
    int grid_y_start = (int)(particle_y / grid_interval);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    Eigen::Vector4d bilinear_coeff;
    get_bilinear_coeff(bilinear_coeff, grid_interval, grid_x_start, grid_x_end, grid_y_start, grid_y_end, particle_x, particle_y);

    int i_idx, j_idx;
    get_matrix_index_2d(grid_x_start, grid_y_start, len_x, len_y+1, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y) && (j_idx >= 0) && (j_idx <= len_x - 1))
        M_v(i_idx, j_idx) += bilinear_coeff[0] * v_particle;

    get_matrix_index_2d(grid_x_start, grid_y_end, len_x, len_y+1, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y) && (j_idx >= 0) && (j_idx <= len_x - 1))
        M_v(i_idx, j_idx) += bilinear_coeff[1] * v_particle;

    get_matrix_index_2d(grid_x_end, grid_y_start, len_x, len_y+1, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y) && (j_idx >= 0) && (j_idx <= len_x - 1))
        M_v(i_idx, j_idx) += bilinear_coeff[2] * v_particle;

    get_matrix_index_2d(grid_x_end, grid_y_end, len_x, len_y+1, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y) && (j_idx >= 0) && (j_idx <= len_x - 1))
        M_v(i_idx, j_idx) += bilinear_coeff[3] * v_particle;
}
