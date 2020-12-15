#include "grid.h"

/*
    Compute gradient of the pressure at given grid location (x,y), storing to gradient_p
    EFFECT: Updates gradient_p
*/
void get_gradient_pressure(Eigen::Vector4d gradient_p, Eigen::MatrixXd pressures, int x, int y, double interval_x, double interval_y, double alpha)
{
    assert(interval_x > 0 && "interval x should be greater than 0.");
    assert(interval_y > 0 && "interval y should be greater than 0.");
    assert(alpha != 1 && "alpha should not be 1.");

    Eigen::Matrix4d D;
    D.setZero();
    double inverse_int_x = 1.0 / interval_x;
    double inverse_int_y = 1.0 / interval_y;
    /* with ghost pressure */
    D(0, 1) = inverse_int_x * (1.0 + (alpha / (1.0 - alpha)));
    D(1, 0) = inverse_int_x;
    D(1, 1) = -inverse_int_x;
    D(2, 1) = inverse_int_y;
    D(2, 2) = -inverse_int_y;
    D(3, 1) = -inverse_int_y;
    D(3, 3) = inverse_int_y;

    Eigen::Vector4d p(pressures(x + 2, y), pressures(x + 1, y), pressures(x + 1, y - 1), pressures(x + 1, y + 1));
    gradient_p = D * p;
}

/*
    Accumulate v of the single particle onto the whole grid matrix using bilinear interpolation.
    There are four nodes affected by the particle, four corners correspondingly.
    interval_x and interval_y are the intervals of grids.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    EFFECT: Updates M_v
*/
void v_particle_onto_grid_v(Eigen::MatrixXd M_v, Eigen::Vector2d pos_particle, double v_particle, double interval_x, double interval_y)
{
    assert(interval_x * interval_y != 0 && "intervals should not be 0.");
    /* pos_particle -> (x,y) */
    double particle_x = pos_particle(0);
    double particle_y = pos_particle(1);
    int grid_x_start = (int)(particle_x / interval_x);
    int grid_y_start = (int)(particle_y / interval_y);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    M_v(grid_x_start, grid_y_start) = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (grid_y_end * interval_y - particle_y) * v_particle;
    M_v(grid_x_start, grid_y_end) = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (particle_y - grid_y_start * interval_y) * v_particle;
    M_v(grid_x_end, grid_y_start) = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (grid_y_end * interval_y - particle_y) * v_particle;
    M_v(grid_x_end, grid_y_end) = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (particle_y - grid_y_start * interval_y) * v_particle;
}

/*
    Accumulate u of the single particle onto the whole grid matrix using bilinear interpolation.
    There are four nodes affected by the particle, four corners correspondingly.
    interval_x and interval_y are the intervals of grids.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    NOTE: the same logic as v_particle_onto_grid_v, can merge into one if we want.
    EFFECT: Updates M_u
*/
void u_particle_onto_grid_u(Eigen::MatrixXd M_u, Eigen::Vector2d pos_particle, double u_particle, double interval_x, double interval_y) {
    assert(interval_x * interval_y != 0 && "intervals should not be 0.");
    /* pos_particle -> (x,y) */
    double particle_x = pos_particle(0);
    double particle_y = pos_particle(1);
    int grid_x_start = (int)(particle_x / interval_x);
    int grid_y_start = (int)(particle_y / interval_y);
    int grid_x_end = grid_x_start + 1;
    int grid_y_end = grid_y_start + 1;

    M_u(grid_x_start, grid_y_start) = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (grid_y_end * interval_y - particle_y) * u_particle;
    M_u(grid_x_start, grid_y_end) = (1 / (interval_x * interval_y)) * (grid_x_end * interval_x - particle_x) * (particle_y - grid_y_start * interval_y) * u_particle;
    M_u(grid_x_end, grid_y_start) = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (grid_y_end * interval_y - particle_y) * u_particle;
    M_u(grid_x_end, grid_y_end) = (1 / (interval_x * interval_y)) * (particle_x - grid_x_start * interval_x) * (particle_y - grid_y_start * interval_y) * u_particle;
}
