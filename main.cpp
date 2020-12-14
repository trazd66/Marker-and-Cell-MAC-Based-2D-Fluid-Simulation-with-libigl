#include <iostream>
#include <thread>

#include <visualization.h>
#include <igl/edges.h>
#include <igl/edge_lengths.h>
#include <igl/readMESH.h>
#include <igl/boundary_facets.h>
#include <Eigen/Dense>

#include <init_state.h>

//Simulation State
bool simulating = true;

double t = 0; //simulation time 
double dt = 0.005; //time step


void simulate() {

    while(simulating) {

    }
}

bool draw(igl::opengl::glfw::Viewer & viewer) {

    return false;
}

int main(int argc, char **argv) {



    //run simulation in seperate thread to avoid slowing down the UI
    std::thread simulation_thread(simulate);
    simulation_thread.detach();

    //setup libigl viewer and activate 
    // Visualize::setup(q, qdot, true);
    // Visualize::add_object_to_scene(V,F, Eigen::RowVector3d(244,165,130)/255.);
    Visualize::viewer().callback_post_draw = &draw;
    
    Visualize::viewer().launch_init(true,false,"fluid-sim",0,0);
    Visualize::viewer().launch_rendering(true);
    simulating = false;
    Visualize::viewer().launch_shut();

}
