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
