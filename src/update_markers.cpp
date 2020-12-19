#include<update_markers.h>
#include<utilities.h>

void update_markers_2d(Eigen::MatrixXd &M_particles,std::vector<int> marker_index, double grid_interval,
                    Eigen::MatrixXd &M_fluid){

                    M_fluid.setZero();
                    for (int i = 0; i < marker_index.size(); i++)
                    {
                        int i_idx,j_idx;
                        get_matrix_index_2d((int)M_particles(marker_index[i], 0)/grid_interval,(int)M_particles(marker_index[i], 1)/grid_interval,M_fluid.cols(),M_fluid.rows(),i_idx,j_idx);        
                        M_fluid(i_idx, j_idx) = 1;
                    }
            
    }