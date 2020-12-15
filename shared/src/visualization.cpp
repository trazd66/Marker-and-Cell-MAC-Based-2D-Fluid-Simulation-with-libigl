#include <visualization.h> 

//libigl viewer
namespace Visualize {

    igl::opengl::glfw::Viewer g_viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    
    //meshes in the scene 
    std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi> > g_geometry;
    std::vector<unsigned int> g_id; //id into libigl for these meshes 

}

igl::opengl::glfw::Viewer & Visualize::viewer() { return g_viewer; }



void Visualize::setup(const Eigen::VectorXd &q, const Eigen::VectorXd &qdot, bool ps_plot) {


    //add new menu for phase space plotting
    Visualize::g_viewer.plugins.push_back(&menu);

    menu.callback_draw_viewer_menu = [&]()
    {
        ImGuiStyle& style = ImGui::GetStyle();
        style.WindowRounding = 5.3f;
        style.FrameRounding = 2.3f;
        style.ScrollbarRounding = 0;

        style.Colors[ImGuiCol_Text]                  = ImVec4(0.0f, 0.0f, 0.0f, 1.0f);
        style.Colors[ImGuiCol_TextDisabled]          = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
        style.Colors[ImGuiCol_WindowBg]              = ImVec4(0.8f, 0.8f, 0.8f, 1.00f);
        style.Colors[ImGuiCol_PopupBg]               = ImVec4(0.05f, 0.05f, 0.10f, 0.85f);
        style.Colors[ImGuiCol_Border]                = ImVec4(0.70f, 0.70f, 0.70f, 0.65f);
        style.Colors[ImGuiCol_BorderShadow]          = ImVec4(1.00f, 0.00f, 0.00f, 0.00f);
        style.Colors[ImGuiCol_FrameBg]               = ImVec4(1.00f, 1.00f, 1.00f, 1.00f);
        style.Colors[ImGuiCol_FrameBgHovered]        = ImVec4(0.90f, 0.80f, 0.80f, 0.40f);
        style.Colors[ImGuiCol_FrameBgActive]         = ImVec4(0.90f, 0.65f, 0.65f, 0.45f);
        style.Colors[ImGuiCol_TitleBg]               = ImVec4(0.960f, 0.960f, 0.960f, 1.0f);
        style.Colors[ImGuiCol_TitleBgCollapsed]      = ImVec4(0.960f, 0.960f, 0.960f, 1.0f);
        style.Colors[ImGuiCol_TitleBgActive]         = ImVec4(0.960f, 0.960f, 0.960f, 1.0f);
        style.Colors[ImGuiCol_MenuBarBg]             = ImVec4(0.01f, 0.01f, 0.02f, 0.80f);
        style.Colors[ImGuiCol_ScrollbarBg]           = ImVec4(0.20f, 0.25f, 0.30f, 0.60f);
        style.Colors[ImGuiCol_ScrollbarGrab]         = ImVec4(0.55f, 0.53f, 0.55f, 0.51f);
        style.Colors[ImGuiCol_ScrollbarGrabHovered]  = ImVec4(0.56f, 0.56f, 0.56f, 1.00f);
        style.Colors[ImGuiCol_ScrollbarGrabActive]   = ImVec4(0.56f, 0.56f, 0.56f, 0.91f);
        style.Colors[ImGuiCol_CheckMark]             = ImVec4(0.90f, 0.90f, 0.90f, 0.83f);
        style.Colors[ImGuiCol_SliderGrab]            = ImVec4(0.70f, 0.70f, 0.70f, 0.62f);
        style.Colors[ImGuiCol_SliderGrabActive]      = ImVec4(0.30f, 0.30f, 0.30f, 0.84f);
        style.Colors[ImGuiCol_Button]                = ImVec4(0.48f, 0.72f, 0.89f, 1.00f);
        style.Colors[ImGuiCol_ButtonHovered]         = ImVec4(0.50f, 0.69f, 0.99f, 1.00f);
        style.Colors[ImGuiCol_ButtonActive]          = ImVec4(0.80f, 0.50f, 0.50f, 1.00f);
        style.Colors[ImGuiCol_Header]                = ImVec4(0.44f, 0.61f, 0.86f, 1.00f);
        style.Colors[ImGuiCol_HeaderHovered]         = ImVec4(0.44f, 0.61f, 0.86f, 1.00f);
        style.Colors[ImGuiCol_HeaderActive]          = ImVec4(0.44f, 0.61f, 0.86f, 1.00f);
        style.Colors[ImGuiCol_ResizeGrip]            = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
        style.Colors[ImGuiCol_ResizeGripHovered]     = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
        style.Colors[ImGuiCol_ResizeGripActive]      = ImVec4(0.70f, 0.70f, 0.70f, 1.00f);
        style.Colors[ImGuiCol_PlotLines]             = ImVec4(0.00f, 1.00f, 0.00f, 1.00f);
        style.Colors[ImGuiCol_PlotLinesHovered]      = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
        style.Colors[ImGuiCol_PlotHistogram]         = ImVec4(0.90f, 0.70f, 0.00f, 1.00f);
        style.Colors[ImGuiCol_PlotHistogramHovered]  = ImVec4(1.00f, 0.60f, 0.00f, 1.00f);
        style.Colors[ImGuiCol_TextSelectedBg]        = ImVec4(0.00f, 0.00f, 1.00f, 0.35f);
        style.Colors[ImGuiCol_ModalWindowDarkening]  = ImVec4(0.20f, 0.20f, 0.20f, 0.35f);

        // Draw parent menu content
        menu.draw_viewer_menu();
    };

    Visualize::g_viewer.core().background_color.setConstant(1.0);
    Visualize::g_viewer.core().is_animating = true;
}

void Visualize::add_object_to_scene(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::RowVector3d color) {

    //add mesh to libigl and store id for access later
    if(g_geometry.size() == 0) {
        g_id.push_back(0);     
    } else {
        g_id.push_back(g_viewer.append_mesh());
    }

    g_viewer.data().set_mesh(V,F);
    g_viewer.data().set_colors(color);

    //add mesh to geometry vector
    g_geometry.push_back(std::make_pair(V,F));

}

void Visualize::rigid_transform_1d(unsigned int id, double x) {
    
    //reset vertex positions 
    for(unsigned int ii=0; ii<g_geometry[id].first.rows(); ++ii) {
        g_viewer.data_list[g_id[id]].V(ii,0) = g_geometry[id].first(ii,0) + x;
        }

        //tell viewer to update
        g_viewer.data_list[g_id[id]].dirty |= igl::opengl::MeshGL::DIRTY_POSITION;
}

void Visualize::scale_x(unsigned int id, double x) {
    
    //reset vertex positions 
    for(unsigned int ii=0; ii<g_geometry[id].first.rows(); ++ii) {
        g_viewer.data_list[g_id[id]].V(ii,0) = x*g_geometry[id].first(ii,0);
        }

        //tell viewer to update
        g_viewer.data_list[g_id[id]].dirty |= igl::opengl::MeshGL::DIRTY_POSITION;

}

void Visualize::update_vertex_positions(unsigned int id, Eigen::Ref<const Eigen::VectorXd> pos) {

    //update vertex positions
    for(unsigned int ii=0; ii<g_geometry[id].first.rows(); ++ii) {
        g_viewer.data_list[g_id[id]].V.row(ii) = pos.segment<3>(3*ii).transpose();
        }

        //tell viewer to update
        g_viewer.data_list[g_id[id]].dirty |= igl::opengl::MeshGL::DIRTY_POSITION;
}


