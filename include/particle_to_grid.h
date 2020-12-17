
/*
    Functions for converting particle verlocities to grid velocities.
*/
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    Accumulate v of the single particle onto the whole grid matrix using bilinear interpolation.
    There are four nodes affected by the particle, four corners correspondingly.
    interval_x and interval_y are the intervals of grids.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    EFFECT: Updates M_v
*/
void v_particle_onto_grid_v(Eigen::MatrixXd &M_v, Eigen::Vector2d &pos_particle, double v_particle, double interval_x, double interval_y, const int len_x,const int len_y);

/*
    Accumulate u of the single particle onto the whole grid matrix using bilinear interpolation.
    There are four nodes affected by the particle, four corners correspondingly.
    interval_x and interval_y are the intervals of grids.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    EFFECT: Updates M_u
    NOTE: the same logic as v_particle_onto_grid_v, can merge into one if we want.
*/
void u_particle_onto_grid_u(Eigen::MatrixXd &M_u, Eigen::Vector2d &pos_particle, double u_particle, double interval_x, double interval_y, const int len_x,const int len_y);
