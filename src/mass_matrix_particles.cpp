#include <mass_matrix_particles.h>

void mass_matrix_particles(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, double mass) {
    M = Eigen::SparseMatrixd();
    M.resize(q.size(),q.size());
    M.setIdentity();
    M *= mass;
}
