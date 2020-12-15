#include <init_state.h>


void init_state_3d(const int bb_size_x,const int bb_size_y,const int bb_size_z, 
                const double grid_interval, const int num_particles, Eigen::MatrixXd &M_particles,
                Eigen::TensorXd &M_u, Eigen::TensorXd &M_v, Eigen::TensorXd &M_w){
                    
                    M_particles = Eigen::MatrixXd(num_particles,3);

                    M_u = Eigen::TensorXd();
                    M_u.init_tensor(bb_size_x + 1,bb_size_y,bb_size_z);

                    M_v = Eigen::TensorXd();
                    M_v.init_tensor(bb_size_x,bb_size_y + 1,bb_size_z);

                    M_w = Eigen::TensorXd();
                    M_w.init_tensor(bb_size_x,bb_size_y,bb_size_z + 1);
                }



void init_state_2d(const int bb_size_x,const int bb_size_y, 
                const double grid_interval, const int num_particles, Eigen::MatrixXd &M_particles,
                Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v,
                Eigen::MatrixXd &M_pressures){
                    // TODO: initialize particle position at random positions within boundaries
                    M_particles = Eigen::MatrixXd(num_particles,2);
                    M_particles_u = Eigen::VectorXd(num_particles);
                    M_particles_v = Eigen::VectorXd(num_particles);

                    M_pressures = Eigen::MatrixXd(bb_size_x,bb_size_y);
                    M_u = Eigen::MatrixXd(bb_size_x + 1,bb_size_y);
                    M_v = Eigen::MatrixXd(bb_size_x,bb_size_y + 1);

                    M_u.setZero();
                    M_v.setZero();
                    M_pressures.setZero();
                    M_particles_u.setZero();
                    M_particles_v.setZero();
                }
