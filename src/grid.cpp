#include "grid.h"

Eigen::MatrixXd InitGridPressures(const int rows, const int cols)
{
    Eigen::MatrixXd result(rows, cols);
    result.setZero();
    return result;
}

Eigen::MatrixXd InitGridVelocityU(const int rows, const int cols)
{
    /* need to add one since num_borders = num_grid + 1 */
    Eigen::MatrixXd result(rows, cols + 1);
    result.setZero();
    return result;
}

Eigen::MatrixXd InitGridVelocityV(const int rows, const int cols)
{
    /* need to add one since num_borders = num_grid + 1 */
    Eigen::MatrixXd result(rows + 1, cols);
    result.setZero();
    return result;
}
