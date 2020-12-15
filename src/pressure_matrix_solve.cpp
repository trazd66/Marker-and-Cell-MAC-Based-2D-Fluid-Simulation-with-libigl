#include <pressure_matrix_solve.h>

/*
    Solves pressure from Ap = f using Conjugate Gradient Method.
    EFFECT: updates pressure p
*/
void solve_pressure_p (Eigen::VectorXd &p, Eigen::MatrixXd A, Eigen::VectorXd f) {
    Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower|Eigen::Upper> cg;
    cg.compute(A);
    p = cg.solve(f);
}
