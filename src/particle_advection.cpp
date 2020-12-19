#include <particle_advection.h>


void advect_particle_2d(Eigen::MatrixXd &M_particles, Eigen::VectorXd &u, Eigen::VectorXd &v, double dt){

    for (int i = 0; i < M_particles.rows(); i++) {
        Eigen::Vector2d particle_pos, particle_velocity;
        particle_velocity << u[i], v[i];
        particle_pos = M_particles.row(i).transpose();
        Eigen::Vector2d k1;
        Eigen::Vector2d k2;
        Eigen::Vector2d k3;
        Eigen::Vector2d k4;
        k1 = particle_pos + dt * particle_velocity / 2;
        k2 = particle_pos + (dt * k1) / 2;
        k3 = particle_pos + (dt * k2) / 2;
        k4 = particle_pos + dt * k3;

        particle_pos += (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        M_particles.row(i) = k1;    
    }
}
