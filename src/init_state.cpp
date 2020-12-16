#include <init_state.h>
#include <stdlib.h>

void init_state_3d(const int bb_size_x, const int bb_size_y, const int bb_size_z,
                   const double grid_interval, const int num_particles, Eigen::MatrixXd &M_particles,
                   Eigen::MatrixXd &M_signed_distance,
                   Eigen::TensorXd &M_u, Eigen::TensorXd &M_v, Eigen::TensorXd &M_w)
{

    M_particles = Eigen::MatrixXd(num_particles, 3);

    M_u = Eigen::TensorXd();
    M_u.init_tensor(bb_size_x + 1, bb_size_y, bb_size_z);

    M_v = Eigen::TensorXd();
    M_v.init_tensor(bb_size_x, bb_size_y + 1, bb_size_z);

    M_w = Eigen::TensorXd();
    M_w.init_tensor(bb_size_x, bb_size_y, bb_size_z + 1);
}

void init_state_2d(const int bb_size_x, const int bb_size_y,
                   const double grid_interval, const int num_particles,
                   Eigen::MatrixXd &M_particles, Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v,
                   Eigen::MatrixXd &M_signed_distance,
                   Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v, Eigen::MatrixXd &M_pressures)
{
    M_particles = Eigen::MatrixXd(num_particles, 2);
    M_particles_u = Eigen::VectorXd(num_particles);
    M_particles_v = Eigen::VectorXd(num_particles);

    M_pressures = Eigen::MatrixXd(bb_size_x, bb_size_y);
    M_u = Eigen::MatrixXd(bb_size_x + 1, bb_size_y);
    M_v = Eigen::MatrixXd(bb_size_x, bb_size_y + 1);

    M_u.setZero();
    M_v.setZero();
    M_pressures.setZero();
    M_particles_u.setZero();
    M_particles_v.setZero();

    /* initialize particle position at random positions within boundaries */
    const int boundary_x = (int)(bb_size_x * grid_interval);
    const int boundary_y = (int)(bb_size_x * grid_interval);
    for (int i = 0; i < num_particles; i++)
    {
        const int boundary_x = (int)(bb_size_x * grid_interval + 1);
        const int boundary_y = (int)(bb_size_y * grid_interval + 1);

        M_particles(i, 0) = (double)(rand() % boundary_x);
        M_particles(i, 1) = (double)(rand() % boundary_y);
    }

    //calculating initial signed distance matrix
    //set to -1 for now
    M_signed_distance = Eigen::MatrixXd(bb_size_x,bb_size_y);
    M_signed_distance.setOnes();
    M_signed_distance *= -1;//everything is inside surface
    //except the boundaries, they should be outside of our surface at all times
    for (int i = 0; i < bb_size_y; i++)
    {
        M_signed_distance(0,i) = 1;
        M_signed_distance(bb_size_x - 1,i) = 1;
    }

    for (int i = 0; i < bb_size_x; i++)
    {
        M_signed_distance(i,0) = 1;
        M_signed_distance(i,bb_size_y - 1) = 1;
    }
        
}
