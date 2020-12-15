/*
    Functions for manipulating grids.
*/
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    Compute gradient of the pressure at given grid location (x,y), storing to gradient_p
    EFFECT: Updates gradient_p
*/
void get_gradient_pressure(Eigen::Vector4d gradient_p, Eigen::MatrixXd pressures, int x, int y, double interval_x, double interval_y, double alpha);

