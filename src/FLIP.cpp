#include <FLIP.h>
#include <utilities.h>
#include <iostream>

/*
    Reconstruct particle velocity u from grid velocity using bilinear interpolation.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    Returns updated u velocity for the particle.
*/
double grid_to_particle_FLIP_u (Eigen::MatrixXd &old_M_u, Eigen::MatrixXd &M_u, Eigen::Vector2d &pos_particle, double particle_u, double grid_interval, const int len_x,const int len_y) {
    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)(particle_x / grid_interval);
    int grid_y_start = (int)((particle_y - 0.5 * grid_interval) / grid_interval);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    Eigen::MatrixXd delta_M_u;
    delta_M_u = M_u - old_M_u;

    // std::cout << "*** PIC_u -> " << grid_x_start << "->" << grid_x_end
    //             << " " << grid_y_start << "->" << grid_y_end
    //             << " len_x: " << len_x << " len_y: " << len_y <<'\n';
    int i_idx, j_idx;
    get_matrix_index_2d(grid_x_start, grid_y_start, len_x+1, len_y, i_idx, j_idx);
    // std::cout << "*** PIC_u -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    double x_start_y_start = (1 / (grid_interval * grid_interval)) * (grid_x_end * grid_interval - particle_x) * (grid_y_end * grid_interval - particle_y) * delta_M_u(i_idx, j_idx);

    get_matrix_index_2d(grid_x_start, grid_y_end, len_x+1, len_y, i_idx, j_idx);
    // std::cout << "*** PIC_u -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    double x_start_y_end = (1 / (grid_interval * grid_interval)) * (grid_x_end * grid_interval - particle_x) * (particle_y - grid_y_start * grid_interval) * delta_M_u(i_idx, j_idx);

    get_matrix_index_2d(grid_x_end, grid_y_start, len_x+1, len_y, i_idx, j_idx);
    // std::cout << "*** PIC_u -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    double x_end_y_start = (1 / (grid_interval * grid_interval)) * (particle_x - grid_x_start * grid_interval) * (grid_y_end * grid_interval - particle_y) * delta_M_u(i_idx, j_idx);

    get_matrix_index_2d(grid_x_end, grid_y_end, len_x+1, len_y, i_idx, j_idx);
    // std::cout << "*** PIC_u -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    double x_end_y_end = (1 / (grid_interval * grid_interval)) * (particle_x - grid_x_start * grid_interval) * (particle_y - grid_y_start * grid_interval) * delta_M_u(i_idx, j_idx);

    double result = particle_u + (x_start_y_start + x_start_y_end + x_end_y_start + x_end_y_end);

    // std::cout << "FLIP u ->" << result << std::endl;
    return result;
}

/*
    Reconstruct particle velocity v from grid velocity using bilinear interpolation.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    Returns updated v velocity for the particle.
*/
double grid_to_particle_FLIP_v (Eigen::MatrixXd &old_M_v, Eigen::MatrixXd &M_v, Eigen::Vector2d &pos_particle, double particle_v, double grid_interval, const int len_x,const int len_y) {
    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)((particle_x - 0.5 * grid_interval) / grid_interval);
    int grid_y_start = (int)(particle_y / grid_interval);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    Eigen::MatrixXd delta_M_v;
    delta_M_v = M_v - old_M_v;

    // std::cout << "*** PIC_v -> " << grid_x_start << "->" << grid_x_end
    //             << " " << grid_y_start << "->" << grid_y_end
    //             << " len_x: " << len_x << " len_y: " << len_y <<'\n';
    int i_idx, j_idx;
    get_matrix_index_2d(grid_x_start, grid_y_start, len_x, len_y+1, i_idx, j_idx);
    // std::cout << "*** PIC_v -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    double x_start_y_start = (1 / (grid_interval * grid_interval)) * (grid_x_end * grid_interval - particle_x) * (grid_y_end * grid_interval - particle_y) * delta_M_v(i_idx, j_idx);

    get_matrix_index_2d(grid_x_start, grid_y_end, len_x, len_y+1, i_idx, j_idx);
    // std::cout << "*** PIC_v -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    double x_start_y_end = (1 / (grid_interval * grid_interval)) * (grid_x_end * grid_interval - particle_x) * (particle_y - grid_y_start * grid_interval) * delta_M_v(i_idx, j_idx);

    get_matrix_index_2d(grid_x_end, grid_y_start, len_x, len_y+1, i_idx, j_idx);
    // std::cout << "*** PIC_v -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    double x_end_y_start = (1 / (grid_interval * grid_interval)) * (particle_x - grid_x_start * grid_interval) * (grid_y_end * grid_interval - particle_y) * delta_M_v(i_idx, j_idx);

    get_matrix_index_2d(grid_x_end, grid_y_end, len_x, len_y+1, i_idx, j_idx);
    // std::cout << "*** PIC_v -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    double x_end_y_end = (1 / (grid_interval * grid_interval)) * (particle_x - grid_x_start * grid_interval) * (particle_y - grid_y_start * grid_interval) * delta_M_v(i_idx, j_idx);

    double result = particle_v + (x_start_y_start + x_start_y_end + x_end_y_start + x_end_y_end);

    // std::cout << "FLIP v ->" << result << std::endl;

    return result;
}
