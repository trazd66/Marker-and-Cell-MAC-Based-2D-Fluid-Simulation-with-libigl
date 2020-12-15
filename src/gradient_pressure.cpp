#include <gradient_pressure.h>

/*
    Compute gradient of the pressure at given grid location (x,y), storing to gradient_p
    EFFECT: Updates gradient_p
*/
void get_gradient_pressure(Eigen::Vector4d &gradient_p, Eigen::MatrixXd pressures, int x, int y, double interval_x, double interval_y, double alpha)
{
    assert(interval_x > 0 && "interval x should be greater than 0.");
    assert(interval_y > 0 && "interval y should be greater than 0.");
    assert(alpha != 1 && "alpha should not be 1.");

    Eigen::Matrix4d D;
    D.setZero();
    double inverse_int_x = 1.0 / interval_x;
    double inverse_int_y = 1.0 / interval_y;
    /* with ghost pressure */
    D(0, 1) = inverse_int_x * (1.0 + (alpha / (1.0 - alpha)));
    D(1, 0) = inverse_int_x;
    D(1, 1) = -inverse_int_x;
    D(2, 1) = inverse_int_y;
    D(2, 2) = -inverse_int_y;
    D(3, 1) = -inverse_int_y;
    D(3, 3) = inverse_int_y;

    Eigen::Vector4d p(pressures(x + 2, y), pressures(x + 1, y), pressures(x + 1, y - 1), pressures(x + 1, y + 1));
    gradient_p = D * p;
}
