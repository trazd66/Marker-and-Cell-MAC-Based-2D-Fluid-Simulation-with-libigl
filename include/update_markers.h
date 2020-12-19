#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>


/***
 * 
 * Given the u and v velocity grid, update the M_liquid matrix from all the marker particles
 * 
 ***/
void update_markers_2d(Eigen::MatrixXd &M_particles,std::vector<int> marker_index, double grid_interval,
                    Eigen::MatrixXd &M_fluid);
