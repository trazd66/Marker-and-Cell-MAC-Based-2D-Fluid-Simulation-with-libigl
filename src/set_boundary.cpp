#include <set_boundary.h>

void set_boundary_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v)
{
    for (int i = 0; i < M_u.rows(); i++)
    {
        M_u(i, 0) = 0;
        M_u(i, M_u.rows()) = 0;
    }
    for (int i = 0; i < M_v.cols(); i++)
    {
        M_v(0, i) = 0;
        M_v(M_v.cols(), i) = 0;
    }
}
