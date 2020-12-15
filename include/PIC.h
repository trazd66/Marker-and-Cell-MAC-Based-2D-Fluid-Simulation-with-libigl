/*
    Particle-In-Cell Transfer for reconstructing particle velocity from grid velocities.
*/

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    Reconstruct particle velocity u from grid velocity using bilinear interpolation.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    Returns updated u velocity for the particle.
*/
double grid_to_particle_PIC_u (Eigen::MatrixXd M_u, Eigen::Vector2d pos_particle, double interval_x, double interval_y);
/*
    Reconstruct particle velocity v from grid velocity using bilinear interpolation.
    https://en.wikipedia.org/wiki/Bilinear_interpolation#:~:text=In%20mathematics%2C%20bilinear%20interpolation%20is,again%20in%20the%20other%20direction.
    Returns updated v velocity for the particle.
*/
double grid_to_particle_PIC_v (Eigen::MatrixXd M_v, Eigen::Vector2d pos_particle, double interval_x, double interval_y);
