#include <pressure_matrix_solve.h>
#include <iostream>

/*
    Solves pressure from Ap = f using Conjugate Gradient Method.
    EFFECT: updates pressure p
*/
void solve_pressure_p (Eigen::VectorXd &p, Eigen::SparseMatrixd &A, Eigen::VectorXd &f) {
    // std::cout << "A empty? -> " << A.nonZeros() << "\n";
    // std::cout << "f empty? -> " << f.isZero() << "\n";

    Eigen::ConjugateGradient<Eigen::SparseMatrixd, Eigen::Lower|Eigen::Upper> cg;
    cg.compute(A);
    p = cg.solve(f);

    // std::cout << "p empty? -> " << p.isZero() << "\n";
}
