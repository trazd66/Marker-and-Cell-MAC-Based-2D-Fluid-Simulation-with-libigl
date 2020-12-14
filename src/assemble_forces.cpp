#include <assemble_forces.h>
#include <dV_spring_particle_particle_dq.h>
#include <mass_matrix_particles.h>
#include <iostream>

void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double mass, double k) { 
        
        f = Eigen::VectorXd(q.size());
        f.setZero();
        Eigen::Vector3d v0,v1;
        Eigen::Vector6d f_i;
        f_i.setZero();
        for (int i = 0; i < E.rows(); i++)
        {
            v0 << q(3*E(i,0)),q(3*E(i,0) + 1),q(3*E(i,0) + 2);
            v1 << q(3*E(i,1)),q(3*E(i,1) + 1),q(3*E(i,1) + 2);
            dV_spring_particle_particle_dq(f_i,v0,v1,l0(i),k);

            f.segment(3*E(i,0),3) += f_i.head(3);
            f.segment(3*E(i,1),3) += f_i.tail(3);
        }
        if (abs(f.sum()) < 5e-06) f.setZero();
        f *= -1;

        std::cout <<"forces: " << f.sum() <<'\n';
        

    };