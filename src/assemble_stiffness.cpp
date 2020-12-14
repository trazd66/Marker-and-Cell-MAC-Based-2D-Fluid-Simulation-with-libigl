#include <assemble_stiffness.h>
#include <d2V_spring_particle_particle_dq2.h>
#include <iostream>
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, 
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> E, Eigen::Ref<const Eigen::VectorXd> l0, 
                     double k) { 
        
        
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(36*E.rows());

        K = Eigen::SparseMatrixd(q.size(),q.size());
        K.setZero();
        Eigen::Vector3d v0,v1;
        Eigen::Matrix66d H_i,H_i_T;        
        for (int i = 0; i < E.rows(); i++)
        {
            v0 << q(3*E(i,0)),q(3*E(i,0) + 1),q(3*E(i,0) + 2);
            v1 << q(3*E(i,1)),q(3*E(i,1) + 1),q(3*E(i,1) + 2);

            d2V_spring_particle_particle_dq2(H_i,v0,v1,l0(i),k);
            H_i_T = H_i.transpose().eval();
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    triplets.push_back(Eigen::Triplet<double>(3 * E(i,0) + j,3 * E(i,0) + k ,H_i(0 + j,0 + k)));
                    triplets.push_back(Eigen::Triplet<double>(3 * E(i,0) + j,3 * E(i,1) + k ,H_i(0 + j,3 + k)));
                    triplets.push_back(Eigen::Triplet<double>(3 * E(i,1) + j,3 * E(i,0) + k ,H_i(3 + j,0 + k)));
                    triplets.push_back(Eigen::Triplet<double>(3 * E(i,1) + j,3 * E(i,1) + k ,H_i(3 + j,3 + k)));
                }
            }
            // std::cout << H_i <<'\n' << '\n';
    
        }

        
        K.setFromTriplets(triplets.begin(), triplets.end());

        K *= -1;


    };