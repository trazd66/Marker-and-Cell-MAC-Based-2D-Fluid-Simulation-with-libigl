/*
    Functions for manipulating grids.
*/
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    computes gradient of the pressure at given grid location (x,y)
*/
Eigen::Vector4d get_gradient_pressure(Eigen::MatrixXd pressures, int x, int y, double interval_x, double interval_y, double alpha);
