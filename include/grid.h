/*
    Functions for manipulating grids.
*/
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    return the grid matrix of the grid
    rows: number of grids in rows
    cols: nubmer of grids in columns
*/
Eigen::MatrixXd InitGridPressures(const int rows, const int cols);
/*
    return the u (horizontal) velocity matrix of the grids
    rows: number of grids in rows
    cols: nubmer of grids in columns
*/
Eigen::MatrixXd InitGridVelocityU(const int rows, const int cols);
/*
    return the v (vertical) velocity matrix of the grids
    rows: number of grids in rows
    cols: nubmer of grids in columns
*/
Eigen::MatrixXd InitGridVelocityV(const int rows, const int cols);

/*
    computes gradient of the pressure at given grid location (x,y)
*/
Eigen::Vector4d get_gradient_pressure(Eigen::MatrixXd pressures, int x, int y, double interval_x, double interval_y, double alpha);
