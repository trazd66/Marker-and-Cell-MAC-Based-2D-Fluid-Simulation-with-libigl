#ifndef  VISUALIZATION_H
#define  VISUALIZATION_H

#define IMGUI_DEFINE_MATH_OPERATORS

#include <igl/unproject.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>

//stl
#include <vector>
#include <array>
#include <deque>

//Eigen
#include <Eigen/Dense>

namespace Visualize {

    //custom phase space plot
    bool plot_phase_space(const char *label, ImVec2 q_bounds, ImVec2 q_dot_bounds, const Eigen::VectorXd &q, const Eigen::VectorXd &q_dot);
    
    void setup(const Eigen::VectorXd &q, const Eigen::VectorXd &qdot, bool ps_plot = false);

    void add_object_to_scene(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::RowVector3d color);
    
    //animate geometry using physics simulation
    void rigid_transform_1d(unsigned int id, double x);

    void scale_x(unsigned int id, double x);

    void update_vertex_positions(unsigned int id, Eigen::Ref<const Eigen::VectorXd> pos);


    igl::opengl::glfw::Viewer & viewer();

}


#endif