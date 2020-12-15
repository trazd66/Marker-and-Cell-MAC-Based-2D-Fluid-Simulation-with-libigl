#include <pressure_matrix_solve.h>

void solve_pressure_p (Eigen::VectorXd &p, Eigen::MatrixXd A, Eigen::VectorXd f) {
    Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower|Eigen::Upper> cg;
    cg.compute(A);
    p = cg.solve(f);
}
