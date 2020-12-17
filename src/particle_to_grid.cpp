#include <particle_to_grid.h>
#include <iostream>
#include <utilities.h>
/*
    Accumulate v of the single particle onto the whole grid matrix using bilinear interpolation.
    There are four nodes affected by the particle, four corners correspondingly.
    interval_x and interval_y are the intervals of grids.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    EFFECT: Updates M_v
*/
void v_particle_onto_grid_v(Eigen::MatrixXd &M_v, Eigen::Vector2d &pos_particle, double v_particle, double interval_x, double interval_y, const int len_x,const int len_y)
{
    assert(interval_x * interval_y != 0 && "intervals should not be 0.");
    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)((particle_x - 0.5 * interval_x) / interval_x);
    int grid_y_start = (int)(particle_y / interval_y);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    std::cout << "*** v -> " << grid_x_start << "->" << grid_x_end
                << " " << grid_y_start << "->" << grid_y_end
                << " len_x: " << len_x << " len_y: " << len_y <<'\n';
    int i_idx, j_idx;
    get_matrix_index_2d(grid_x_start, grid_y_start, len_x, len_y+1, i_idx, j_idx);
    std::cout << "*** v -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    M_v(i_idx, j_idx) += (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (grid_y_end * interval_y - particle_y) * v_particle;

    get_matrix_index_2d(grid_x_start, grid_y_end, len_x, len_y+1, i_idx, j_idx);
    std::cout << "*** v -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    M_v(i_idx, j_idx) += (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (particle_y - grid_y_start * interval_y) * v_particle;

    get_matrix_index_2d(grid_x_end, grid_y_start, len_x, len_y+1, i_idx, j_idx);
    std::cout << "*** v -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    M_v(i_idx, j_idx) += (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (grid_y_end * interval_y - particle_y) * v_particle;

    get_matrix_index_2d(grid_x_end, grid_y_end, len_x, len_y+1, i_idx, j_idx);
    std::cout << "*** v -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    M_v(i_idx, j_idx) += (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (particle_y - grid_y_start * interval_y) * v_particle;

}

/*
    Accumulate u of the single particle onto the whole grid matrix using bilinear interpolation.
    There are four nodes affected by the particle, four corners correspondingly.
    interval_x and interval_y are the intervals of grids.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    NOTE: the same logic as v_particle_onto_grid_v, can merge into one if we want.
    EFFECT: Updates M_u
*/
void u_particle_onto_grid_u(Eigen::MatrixXd &M_u, Eigen::Vector2d &pos_particle, double u_particle, double interval_x, double interval_y, const int len_x,const int len_y) {
    assert(interval_x * interval_y != 0 && "intervals should not be 0.");
    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)(particle_x / interval_x);
    int grid_y_start = (int)((particle_y - 0.5 * interval_y) / interval_y);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    std::cout << "*** u -> " << grid_x_start << "->" << grid_x_end
                << " " << grid_y_start << "->" << grid_y_end
                << " len_x: " << len_x << " len_y: " << len_y <<'\n';
    int i_idx, j_idx;
    get_matrix_index_2d(grid_x_start, grid_y_start, len_x+1, len_y, i_idx, j_idx);
    std::cout << "*** u -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    M_u(i_idx, j_idx) += (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (grid_y_end * interval_y - particle_y) * u_particle;

    get_matrix_index_2d(grid_x_start, grid_y_end, len_x+1, len_y, i_idx, j_idx);
    std::cout << "*** u -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    M_u(i_idx, j_idx) += (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (particle_y - grid_y_start * interval_y) * u_particle;

    get_matrix_index_2d(grid_x_end, grid_y_start, len_x+1, len_y, i_idx, j_idx);
    std::cout << "*** u -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    M_u(i_idx, j_idx) += (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (grid_y_end * interval_y - particle_y) * u_particle;

    get_matrix_index_2d(grid_x_end, grid_y_end, len_x+1, len_y, i_idx, j_idx);
    std::cout << "*** u -> " << "i_idx: " << i_idx << " j_idx: " << j_idx << std::endl;
    M_u(i_idx, j_idx) += (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (particle_y - grid_y_start * interval_y) * u_particle;
}
