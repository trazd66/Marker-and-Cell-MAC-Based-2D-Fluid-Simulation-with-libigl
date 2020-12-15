/*
    Functions for manipulating grids.
*/
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    Compute gradient of the pressure at given grid location (x,y), storing to gradient_p
*/
void get_gradient_pressure(Eigen::Vector4d gradient_p, Eigen::MatrixXd pressures, int x, int y, double interval_x, double interval_y, double alpha);

/*
    Accumulate v of the single particle onto the whole grid matrix using bilinear interpolation.
    There are four nodes affected by the particle, four corners correspondingly.
    interval_x and interval_y are the intervals of grids.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
*/
void v_particle_onto_grid_v(Eigen::MatrixXd M_v, Eigen::Vector2d pos_particle, double v_particle, double interval_x, double interval_y);
