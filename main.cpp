#include <iostream>
#include <thread>
#include <string>
#include <unistd.h>

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
#include <utilities.h>
#include <FLIP.h>
#include <external_force_apply.h>
#include <particle_advection.h>
#include <update_markers.h>
#include <velocity_extrapolation.h>
#include <set_boundary.h>

//Simulation State
bool simulating = true;
bool is2d = true;

double rho = 1000.0;    // density of water
double t = 0;           //simulation time
double dt = 0.005;      //time step
double FLIP_potion = 0; // percentage of FLIP result in particle velocity outof (FLIP + PIC)

const int bb_size_x = 50;         // x dimension of the bounding box -> number of grids in x axis
const int bb_size_y = 50;         // y dimension of the bounding box -> number of grids in y axis
const double grid_interval = 0.1; // size of the grid interval, determines the number of grid cells
int num_particles;                // number of particles
Eigen::MatrixXd M_particles;      // the particle matrix
Eigen::VectorXd M_particles_u;    // particle velocity u
Eigen::VectorXd M_particles_v;    // particle velocity v
Eigen::MatrixXd M_u;              // M_u a 2D matrix that contains the x velociies of the grid
Eigen::MatrixXd M_v;              // M_v a 2D matrix that contains the y velociies of the grid
Eigen::MatrixXd M_pressure;       // Grid pressure matrix
Eigen::MatrixXd M_fluid;          // markers for fluid grids, 1 for fluid, 0 for air

void simulate()
{
    // stores grid velocities of last iteration, used by FLIP to calculate delta_u and delta_v
    Eigen::MatrixXd old_M_u;
    Eigen::MatrixXd old_M_v;
    old_M_u = M_u;
    old_M_v = M_v;

    while (simulating)
    {

        /*
            1. Convert particle velocity to grid velocity with bilinear interpolation.
            2. Calculate advection with a second order Runge Kutta method.
            3. Velocity extrapolation from fluid cells onto air cells.
            4. Apply external forces, gravity in our simplified version.
            5. Calculate new pressure using equation from lecture slide: A * p = f
            6. Convert grid velocity back to particle velocity using bilinear interpolation.
        */

        M_u.setZero();
        M_v.setZero();

        // 1.
        for (int i = 0; i < num_particles; i++)
        {
            Eigen::Vector2d particle_pos;
            particle_pos = M_particles.row(i).transpose();
            double u_particle = M_particles_u[i];
            double v_particle = M_particles_v[i];
            v_particle_onto_grid_v(M_v, particle_pos, v_particle, grid_interval, bb_size_x, bb_size_y);
            u_particle_onto_grid_u(M_u, particle_pos, u_particle, grid_interval, bb_size_x, bb_size_y);
        }

        // make sure dt * particle_velocity < grid_interval
        double dt_new = dt;
        update_dt(dt_new, num_particles, grid_interval, M_particles_u, M_particles_v);
    	double num_substeps = std::ceil(dt/dt_new);

        for( int s = 0; s < num_substeps ; s++ ){
            // 2.
            advect_particle_2d(M_particles, M_particles_u, M_particles_v, dt/num_substeps, M_u, M_v, grid_interval, bb_size_x, bb_size_y);
        }

        // clip the particle positions
        // M_particles = M_particles.cwiseMin(grid_interval * bb_size_x).cwiseMax(0);


        // marks fluid cells and air cells
        update_markers_2d(M_particles, grid_interval, M_fluid);

        // 3. extrapolate fluid cell velocity onto air cells
        extrapolate_velocity_2d(M_u, M_v, M_fluid);

        // normalize grid to make sure we have a sane grid velocity
        normalize_grid(M_u, M_v);

        // 4. apply external forces
        apply_external_force_2d(M_u, M_v, dt);

        // make sure boundary cells have zero grid velocities
        set_boundary_2d(M_u, M_v);

        // 5. new pressure calculation using A * p = f
        Eigen::SparseMatrixd A;
        Eigen::VectorXd f;
        assemble_pressure_A_2d(M_u, M_v, M_particles, M_fluid, A);
        assemble_pressure_f_2d(rho, grid_interval, grid_interval, dt, M_u, M_v, M_fluid, M_particles, f);
        grid_pressure_gradient_update_2d(M_u, M_v, M_particles, M_pressure, M_fluid, A, f, rho, dt, grid_interval);

        // 6.
        for (int i = 0; i < num_particles; i++)
        {
            Eigen::Vector2d particle_pos;

            particle_pos = M_particles.row(i).transpose();
            double u_particle = M_particles_u[i];
            double v_particle = M_particles_v[i];

            double new_u_PIC = grid_to_particle_PIC_u(M_u, particle_pos, grid_interval, bb_size_x, bb_size_y);
            double new_v_PIC = grid_to_particle_PIC_v(M_v, particle_pos, grid_interval, bb_size_x, bb_size_y);

            double new_u_FLIP = grid_to_particle_FLIP_u(old_M_u, M_u, particle_pos, u_particle, grid_interval, bb_size_x, bb_size_y);
            double new_v_FLIP = grid_to_particle_FLIP_v(old_M_v, M_v, particle_pos, v_particle, grid_interval, bb_size_x, bb_size_y);

            double new_u = FLIP_potion * new_u_FLIP + (1 - FLIP_potion) * new_u_PIC;
            double new_v = FLIP_potion * new_v_FLIP + (1 - FLIP_potion) * new_v_PIC;

            M_particles_u[i] = new_u;
            M_particles_v[i] = new_v;
        }

        old_M_u = M_u;
        old_M_v = M_v;

        // normalize particle velocity to make sure we have a sane particle velocity
        // normalize_velocity(M_particles_u, M_particles_v);

        t += dt;
    }
}

bool draw(igl::opengl::glfw::Viewer &viewer)
{
    /* water particles */
    Eigen::MatrixXd particle_colors(num_particles, 3);
    particle_colors.setOnes();
    viewer.data().set_points(M_particles, particle_colors);

    /* draw boundaries */
    Eigen::MatrixXd points(4, 2);
    points.row(0) << 0, 0;
    points.row(1) << 0, bb_size_x * grid_interval;
    points.row(2) << bb_size_x * grid_interval, 0;
    points.row(3) << bb_size_x * grid_interval, bb_size_x * grid_interval;

    viewer.data().add_points(points, Eigen::RowVector3d(255, 0, 0));
    viewer.data().add_edges(points.row(0), points.row(1), Eigen::RowVector3d(255, 0, 0));
    viewer.data().add_edges(points.row(1), points.row(3), Eigen::RowVector3d(255, 0, 0));
    viewer.data().add_edges(points.row(3), points.row(2), Eigen::RowVector3d(255, 0, 0));
    viewer.data().add_edges(points.row(2), points.row(0), Eigen::RowVector3d(255, 0, 0));

    /* draw corners */
    std::stringstream l1;
    std::stringstream l2;
    std::stringstream l3;
    std::stringstream l4;
    l1 << "(" << points(0, 0) << ", " << points(0, 1) << ")";
    l1 << "(" << points(1, 0) << ", " << points(1, 1) << ")";
    l1 << "(" << points(2, 0) << ", " << points(2, 1) << ")";
    l1 << "(" << points(3, 0) << ", " << points(3, 1) << ")";
    viewer.data().add_label(points.row(0), l1.str());
    viewer.data().add_label(points.row(1), l2.str());
    viewer.data().add_label(points.row(2), l3.str());
    viewer.data().add_label(points.row(3), l4.str());

    viewer.data().point_size = 10;

    return false;
}

int main(int argc, char **argv)
{
    // initial setup
    init_state_2d(bb_size_x, bb_size_y,
                  grid_interval, num_particles,
                  M_particles,
                  M_fluid,
                  M_u, M_v,
                  M_particles_u, M_particles_v,
                  M_pressure);

    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulate);
    simulation_thread.detach();

    Visualize::viewer().callback_post_draw = &draw;
    Visualize::viewer().launch_init(true, false, "fluid-sim", 0, 0);
    Visualize::viewer().launch_rendering(true);
    simulating = false;
    Visualize::viewer().launch_shut();
}
