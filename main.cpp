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
#include <gradient_pressure.h>
#include <pressure_matrix_solve.h>
#include <PIC.h>

//Simulation State
bool simulating = true;

double t = 0;      //simulation time
double dt = 0.005; //time step

void simulate()
{

    while (simulating)
    {
        /*
            1. Advection natively satisfied because no acceleration involved during movement of particles.
            2. TODO: update particle velocities due to gravity.
            3. TODO: for each particle update grid velocity (particle -> grid).
            4. TODO: for each particle do pressure projection.
            5. TODO: for each particle update particle velocity from grid (grid -> particle).
        */

        t += dt;
    }
}

bool draw(igl::opengl::glfw::Viewer &viewer)
{

    return false;
}

int main(int argc, char **argv)
{

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
