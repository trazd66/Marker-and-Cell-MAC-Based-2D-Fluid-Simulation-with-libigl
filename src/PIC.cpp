#include <PIC.h>

/*
    Reconstruct particle velocity u from grid velocity using bilinear interpolation.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    Returns updated u velocity for the particle.
*/
double grid_to_particle_PIC_u (Eigen::MatrixXd &M_u, Eigen::Vector2d &pos_particle, double interval_x, double interval_y) {
    assert(interval_x * interval_y != 0 && "intervals should not be 0.");
    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)(particle_x / interval_x);
    int grid_y_start = (int)(particle_y / interval_y);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    double x_start_y_start = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (grid_y_end * interval_y - particle_y) * M_u(grid_y_start, grid_x_start);
    double x_start_y_end = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (particle_y - grid_y_start * interval_y) * M_u(grid_y_start, grid_x_end);
    double x_end_y_start = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (grid_y_end * interval_y - particle_y) * M_u(grid_y_end, grid_x_start);
    double x_end_y_end = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (particle_y - grid_y_start * interval_y) * M_u(grid_y_end, grid_x_end);
    double result = x_start_y_start + x_start_y_end + x_end_y_start + x_end_y_end;

    return result;
}

/*
    Reconstruct particle velocity v from grid velocity using bilinear interpolation.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    Returns updated v velocity for the particle.
*/
double grid_to_particle_PIC_v (Eigen::MatrixXd &M_v, Eigen::Vector2d &pos_particle, double interval_x, double interval_y) {
    assert(interval_x * interval_y != 0 && "intervals should not be 0.");
    /* pos_particle -> (x,y) */
    double particle_x = pos_particle[0];
    double particle_y = pos_particle[1];
    int grid_x_start = (int)(particle_x / interval_x);
    int grid_y_start = (int)(particle_y / interval_y);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    double x_start_y_start = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (grid_y_end * interval_y - particle_y) * M_v(grid_y_start, grid_x_start);
    double x_start_y_end = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (particle_y - grid_y_start * interval_y) * M_v(grid_y_start, grid_x_end);
    double x_end_y_start = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (grid_y_end * interval_y - particle_y) * M_v(grid_y_end, grid_x_start);
    double x_end_y_end = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (particle_y - grid_y_start * interval_y) * M_v(grid_y_end, grid_x_end);
    double result = x_start_y_start + x_start_y_end + x_end_y_start + x_end_y_end;

    return result;
}
