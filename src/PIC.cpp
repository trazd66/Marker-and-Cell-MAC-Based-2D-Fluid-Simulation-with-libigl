#include <PIC.h>
#include <utilities.h>

/*
    Reconstruct particle velocity u from grid velocity using bilinear interpolation.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    Returns updated u velocity for the particle.
*/
double grid_to_particle_PIC_u (Eigen::MatrixXd &M_u, Eigen::Vector2d &pos_particle, double interval_x, double interval_y, const int len_x,const int len_y) {

    if (is_out_of_boundary(pos_particle, interval_x, len_x, len_y)) return 0;

    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)(particle_x / interval_x);
    int grid_y_start = (int)((particle_y - 0.5 * interval_y) / interval_y);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    int i_idx, j_idx;
    double x_start_y_start = 0.0;
    double x_start_y_end = 0.0;
    double x_end_y_start = 0.0;
    double x_end_y_end = 0.0;

    get_matrix_index_2d(grid_x_start, grid_y_start, len_x+1, len_y, i_idx, j_idx);
<<<<<<< HEAD
    ASSERT((i_idx >= 0) && (i_idx <= len_y - 1), "i_idx -> " << i_idx << particle_x << ' ' << particle_y);
    ASSERT((j_idx >= 0) && (j_idx <= len_x), "j_idx -> " << j_idx);
    double x_start_y_start = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (grid_y_end * interval_y - particle_y) * M_u(i_idx, j_idx);
=======
    if((i_idx >= 0) && (i_idx <= len_y - 1) && (j_idx >= 0) && (j_idx <= len_x))
        x_start_y_start = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (grid_y_end * interval_y - particle_y) * M_u(i_idx, j_idx);
>>>>>>> c89670e34b0479c7c09fe532b60af5e28a55e63f

    get_matrix_index_2d(grid_x_start, grid_y_end, len_x+1, len_y, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y - 1) && (j_idx >= 0) && (j_idx <= len_x))
        x_start_y_end = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (particle_y - grid_y_start * interval_y) * M_u(i_idx, j_idx);

    get_matrix_index_2d(grid_x_end, grid_y_start, len_x+1, len_y, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y - 1) && (j_idx >= 0) && (j_idx <= len_x))
        x_end_y_start = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (grid_y_end * interval_y - particle_y) * M_u(i_idx, j_idx);

    get_matrix_index_2d(grid_x_end, grid_y_end, len_x+1, len_y, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y - 1) && (j_idx >= 0) && (j_idx <= len_x))
        x_end_y_end = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (particle_y - grid_y_start * interval_y) * M_u(i_idx, j_idx);

    double result = x_start_y_start + x_start_y_end + x_end_y_start + x_end_y_end;

    return result;
}

/*
    Reconstruct particle velocity v from grid velocity using bilinear interpolation.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    Returns updated v velocity for the particle.
*/
double grid_to_particle_PIC_v (Eigen::MatrixXd &M_v, Eigen::Vector2d &pos_particle, double interval_x, double interval_y, const int len_x,const int len_y) {

    if (is_out_of_boundary(pos_particle, interval_x, len_x, len_y)) return 0;

    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)((particle_x - 0.5 * interval_x) / interval_x);
    int grid_y_start = (int)(particle_y / interval_y);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    int i_idx, j_idx;
    double x_start_y_start = 0.0;
    double x_start_y_end = 0.0;
    double x_end_y_start = 0.0;
    double x_end_y_end = 0.0;

    get_matrix_index_2d(grid_x_start, grid_y_start, len_x, len_y+1, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y) && (j_idx >= 0) && (j_idx <= len_x - 1))
        x_start_y_start = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (grid_y_end * interval_y - particle_y) * M_v(i_idx, j_idx);

    get_matrix_index_2d(grid_x_start, grid_y_end, len_x, len_y+1, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y) && (j_idx >= 0) && (j_idx <= len_x - 1))
        x_start_y_end = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (particle_y - grid_y_start * interval_y) * M_v(i_idx, j_idx);

    get_matrix_index_2d(grid_x_end, grid_y_start, len_x, len_y+1, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y) && (j_idx >= 0) && (j_idx <= len_x - 1))
        x_end_y_start = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (grid_y_end * interval_y - particle_y) * M_v(i_idx, j_idx);

    get_matrix_index_2d(grid_x_end, grid_y_end, len_x, len_y+1, i_idx, j_idx);
    if((i_idx >= 0) && (i_idx <= len_y) && (j_idx >= 0) && (j_idx <= len_x - 1))
        x_end_y_end = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (particle_y - grid_y_start * interval_y) * M_v(i_idx, j_idx);

    double result = x_start_y_start + x_start_y_end + x_end_y_start + x_end_y_end;

    return result;
}
