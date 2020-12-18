#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>


/***
 * 
 * Given the u and v velocity grid, update the M_liquid matrix from all the marker particles
 * 
 ***/
void update_markers_2d(Eigen::MatrixXd &M_u,Eigen::MatrixXd &M_v,
                    Eigen::MatrixXd &M_particles,std::vector<int> marker_index,
                    Eigen::MatrixXd &M_liquid){

}