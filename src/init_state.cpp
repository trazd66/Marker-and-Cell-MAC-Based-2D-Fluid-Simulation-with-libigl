#include <init_state.h>
#include <stdlib.h>
#include <utilities.h>
#include <iostream>

void init_state_3d(const int bb_size_x, const int bb_size_y, const int bb_size_z,
                   const double grid_interval, const int num_particles, Eigen::MatrixXd &M_particles,
                   Eigen::MatrixXd &M_fluid,
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

void init_state_2d(const int bb_size_x,const int bb_size_y, 
                const double grid_interval, int &num_particles, 
                Eigen::MatrixXd &M_particles,
                Eigen::MatrixXd &M_fluid,
                Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, 
                Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v, 
                Eigen::MatrixXd &M_pressures)
{

    num_particles = 50 * 20 * 4;
    M_particles = Eigen::MatrixXd(num_particles, 2);
    M_particles_u = Eigen::VectorXd(num_particles);
    M_particles_v = Eigen::VectorXd(num_particles);

    M_pressures = Eigen::MatrixXd(bb_size_y,bb_size_x);
    M_u = Eigen::MatrixXd(bb_size_y,bb_size_x + 1);
    M_v = Eigen::MatrixXd(bb_size_y + 1,bb_size_x);
    M_fluid = Eigen::MatrixXd(bb_size_y,bb_size_x);
    M_u.setZero();
    M_v.setZero();
    M_pressures.setZero();
    M_particles_u.setZero();
    M_particles_v.setZero();
    M_fluid.setOnes(); // initialize all grids as fluid cells

    /* initialize particle position at random positions within boundaries */
    const int boundary_x = (int)(bb_size_x * grid_interval);
    const int boundary_y = (int)(bb_size_y * grid_interval);

    int n = 0;
    for (int x = 0; x < bb_size_x; x++)
    {
        for (int y = 0; y < 20; y++)
        {
            for (int i = 0; i < 4; i++,n++)
            {
                float random_float_x = (rand()%100 / 100.);
                float random_float_y = (rand()%100 / 100.);
                M_particles(n, 0) = (x + random_float_x) * grid_interval;
                M_particles(n, 1) = (y + random_float_y) * grid_interval;
            }
            
        }
    }
    // std::cout << M_particles ;
    // for (int x = 30; x < 5; x++)
    // {
    //     for (int y = 30; y < 5; y++)
    //     {
    //         for (int i = 0; i < 4; i++,n++)
    //         {
    //             float random_float_x = (rand()%100 / 100.);
    //             float random_float_y = (rand()%100 / 100.);
    //             M_particles(n, 0) = (x + random_float_x) * grid_interval;
    //             M_particles(n, 1) = (y + random_float_y) * grid_interval;
    //         }
            
    //     }
    // }

    // for (int i = 0; i < num_particles; i++)
    // {
    //     // const int boundary_x = (int)(bb_size_x * grid_interval + 1);
    //     // const int boundary_y = (int)(bb_size_y * grid_interval + 1);

    //     /*
    //         randome float 0 ~ 1 inclusive, here set to 0.4 ~ 0.6 inclusive
    //         Ref: https://www.softwaretestinghelp.com/random-number-generator-cpp/#:~:text=C%2B%2B%20Random%20Float,-The%20rand%20()&text=We%20can%20generate%20float%20random,and%201.0%20(both%20inclusive).
    //     */
    //     float random_float_x = 0.9 * ((float)rand()/RAND_MAX);
    //     float random_float_y = 0.5 * ((float)rand()/RAND_MAX);
    //     M_particles(i, 0) = (random_float_x * boundary_x);
    //     M_particles(i, 1) = (random_float_y * boundary_y);
    // }
}
