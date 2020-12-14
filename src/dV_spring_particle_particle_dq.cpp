#include <dV_spring_particle_particle_dq.h>
#include <V_spring_particle_particle.h>
void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {
    
    f = Eigen::Vector6d();
    f.setZero();
    double q_1,q_2,q_3,q_4,q_5,q_6;
    q_1 = q0[0];
    q_2 = q0[1];
    q_3 = q0[2];
    q_4 = q1[0];
    q_5 = q1[1];
    q_6 = q1[2];

    f[0] = stiffness*(q_1*2.0-q_4*2.0)*(l0-sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6)))*1.0/sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6))*(-1.0/2.0);
    f[1] = stiffness*(q_2*2.0-q_5*2.0)*(l0-sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6)))*1.0/sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6))*(-1.0/2.0);
    f[2] = stiffness*(q_3*2.0-q_6*2.0)*(l0-sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6)))*1.0/sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6))*(-1.0/2.0);
    f[3] = (stiffness*(q_1*2.0-q_4*2.0)*(l0-sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6)))*1.0/sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6)))/2.0;
    f[4] = (stiffness*(q_2*2.0-q_5*2.0)*(l0-sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6)))*1.0/sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6)))/2.0;
    f[5] = (stiffness*(q_3*2.0-q_6*2.0)*(l0-sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6)))*1.0/sqrt(q_1*(q_1-q_4)+q_2*(q_2-q_5)-q_4*(q_1-q_4)+q_3*(q_3-q_6)-q_5*(q_2-q_5)-q_6*(q_3-q_6)))/2.0;

}