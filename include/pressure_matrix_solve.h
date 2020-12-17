/*
    Ap = f conjugate gradient method to solver for p
*/

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/*
    Solves pressure from Ap = f using Conjugate Gradient Method.
    EFFECT: updates pressure p
*/
void solve_pressure_p (Eigen::VectorXd &p, Eigen::SparseMatrixd &A, Eigen::VectorXd &f);
