#include <V_spring_particle_particle.h>

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {
    Eigen::VectorXd q(6);
    q << q0,q1;
    Eigen::MatrixXd B(3,6);
    B << -1* Eigen::MatrixXd::Identity(3,3),
        Eigen::MatrixXd::Identity(3,3);
    V = 0.5 * stiffness * pow(sqrt(q.transpose() * B.transpose() * B * q) - l0,2);
    
}