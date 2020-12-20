#include <particle_advection.h>
#include <iostream>
#include <PIC.h>
#include <utilities.h>

void advect_particle_2d(Eigen::MatrixXd &M_particles, Eigen::VectorXd &u, Eigen::VectorXd &v, double dt, Eigen::MatrixXd M_u, Eigen::MatrixXd M_v, double grid_interval, int len_x, int len_y){
    for (int i = 0; i < M_particles.rows(); i++) {
        Eigen::Vector2d particle_pos, particle_velocity;
        particle_velocity << u[i], v[i];
        particle_pos = M_particles.row(i).transpose();

        Eigen::Vector2d k1, k2, k3, k4;
        k1 = particle_velocity;

        Eigen::Vector2d k2_velocity;
        Eigen::Vector2d particle_pos_k2 =  particle_pos + 0.5 * dt * k1;
        if (is_out_of_boundary(particle_pos_k2, grid_interval, len_x, len_y)) {
            k2.setZero();
        } else {
            k2[0] = grid_to_particle_PIC_u(M_u, particle_pos_k2, grid_interval, grid_interval, len_x, len_y);
            k2[1] = grid_to_particle_PIC_v(M_v, particle_pos_k2, grid_interval, grid_interval, len_x, len_y);
        }

        Eigen::Vector2d k3_velocity;
        Eigen::Vector2d particle_pos_k3 =  particle_pos + 0.75 * dt * k2;
        if (is_out_of_boundary(particle_pos_k3, grid_interval, len_x, len_y)) {
            k3.setZero();
        } else {
            k3[0] = grid_to_particle_PIC_u(M_u, particle_pos_k3, grid_interval, grid_interval, len_x, len_y);
            k3[1] = grid_to_particle_PIC_v(M_v, particle_pos_k3, grid_interval, grid_interval, len_x, len_y);
        }

        particle_pos += (2.0 / 9.0) * dt * k1 + (3.0 / 9.0) * dt * k2 + (4.0 / 9.0) * dt * k3;
        M_particles.row(i) = particle_pos;
    }
}
