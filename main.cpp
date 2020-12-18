#include <iostream>
#include <thread>

#include <visualization.h>
#include <igl/edges.h>
#include <igl/edge_lengths.h>
#include <igl/readMESH.h>
#include <igl/boundary_facets.h>
#include <Eigen/Dense>

#include <init_state.h>
#include <particle_to_grid.h>
#include <pressure_matrix_solve.h>
#include <PIC.h>
#include <assemble_pressure_A.h>
#include <assemble_pressure_f.h>
#include <grid_pressure_gradient_update.h>

int iteration_counter = 0;

//Simulation State
bool simulating = true;
bool is2d = true;

double rho = 1.0; // density of water
double t = 0;      //simulation time
double dt = 0.005; //time step
Eigen::Vector2d g(0., -0.0098); // gravity acceleration

const int bb_size_x = 50; // x dimension of the bounding box -> number of grids in x axis
const int bb_size_y = 50; // y dimension of the bounding box -> number of grids in y axis
const double grid_interval = 5.; // size of the grid interval, determines the number of grid cells
const int num_particles = 1000; // number of particles
Eigen::MatrixXd M_particles; // the particle matrix
Eigen::VectorXd M_particles_u; // particle velocity u
Eigen::VectorXd M_particles_v; // particle velocity v
Eigen::MatrixXd M_u; // M_u a 2D matrix that contains the x velociies of the grid
Eigen::MatrixXd M_v; // M_v a 2D matrix that contains the y velociies of the grid
Eigen::MatrixXd M_pressures; // Grid pressure matrix
Eigen::MatrixXd M_signed_distance; //Signed_distance matrix

void simulate()
{

    while (simulating)
    {
        std::cout << "iteration " << iteration_counter++ << std::endl;

        /*
            1. Advection natively satisfied because no acceleration involved during movement of particles.
            2. update particle velocities due to gravity.
            3. for each particle update grid velocity (particle -> grid).
            4. for each particle do pressure projection.
            5. for each particle update particle velocity from grid (grid -> particle).
        */

        // 2.
        Eigen::VectorXd g_acc_vector(num_particles);
        g_acc_vector.setOnes();
        M_particles_v += g_acc_vector * g[1] * dt;
        for (int i = 0; i < num_particles; i++) {
           // 3.
           Eigen::Vector2d particle_pos;

           particle_pos = M_particles.row(i).transpose();
           double u_particle = M_particles_u[i];
           double v_particle = M_particles_v[i];
           v_particle_onto_grid_v(M_v, particle_pos, v_particle, grid_interval, grid_interval, bb_size_x, bb_size_y);
           u_particle_onto_grid_u(M_u, particle_pos, u_particle, grid_interval, grid_interval, bb_size_x, bb_size_y);
        }
        std::cout <<"particle_grid_complete"<<'\n';
        // 4.
        Eigen::SparseMatrixd A;
        Eigen::VectorXd f;
        assemble_pressure_A_2d(M_u, M_v, M_particles, M_signed_distance, A);
        std::cout <<"A_complete"<<'\n';

        assemble_pressure_f_2d(rho, grid_interval, grid_interval, dt, M_u, M_v, M_signed_distance, M_particles, f);
        std::cout <<"f_complete"<<'\n';
        grid_pressure_gradient_update_2d(M_u, M_v, M_particles, M_signed_distance, A, f, rho, dt, grid_interval);
        std::cout <<"pressure_complete"<<'\n';
        // 5.
        for (int i = 0; i < num_particles; i++) {
            Eigen::Vector2d particle_pos;

            particle_pos = M_particles.row(i).transpose();
            double u_particle = M_particles_u[i];
            double v_particle = M_particles_v[i];

            double new_u = grid_to_particle_PIC_u (M_u, particle_pos, grid_interval, grid_interval, bb_size_x, bb_size_y);
            double new_v = grid_to_particle_PIC_v (M_v, particle_pos, grid_interval, grid_interval, bb_size_x, bb_size_y);

            M_particles_u[i] = new_u;
            M_particles_v[i] = new_v;
            M_particles(i, 0) += new_u *dt;
            M_particles(i, 1) += new_v *dt;
        }

        t += dt;
    }
}

bool draw(igl::opengl::glfw::Viewer &viewer)
{
    Eigen::MatrixXd particle_colors(num_particles, 3);
    particle_colors.setOnes();
    viewer.data().set_points(M_particles, particle_colors);
    return false;
}

int main(int argc, char **argv)
{
    // initial setup
    init_state_2d(bb_size_x,bb_size_y, 
                grid_interval,num_particles, 
                M_particles,
                M_signed_distance,
                M_u, M_v, 
                M_particles_u, M_particles_v, 
                M_pressures);

    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulate);
    simulation_thread.detach();

    //setup libigl viewer and activate
    // Visualize::setup(q, qdot, true);
    // Visualize::add_object_to_scene(V,F, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::viewer().callback_post_draw = &draw;

    Visualize::viewer().launch_init(true, false, "fluid-sim", 0, 0);
    Visualize::viewer().launch_rendering(true);
    simulating = false;
    Visualize::viewer().launch_shut();
}
