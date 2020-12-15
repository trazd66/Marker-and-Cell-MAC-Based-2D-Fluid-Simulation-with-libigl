#include <init_state.h>


void init_state_3d(const int bb_size_x,const int bb_size_y,const int bb_size_z, 
                const double grid_interval, const int num_particles, Eigen::MatrixXd &M_particles,
                Eigen::TensorXd &M_p, Eigen::TensorXd &M_u, Eigen::TensorXd &M_v, Eigen::TensorXd &M_w){
                    M_particles = Eigen::MatrixXd(num_particles,3);

                    M_p = Eigen::TensorXd();
                    M_p.init_tensor(bb_size_x,bb_size_y,bb_size_z);


                    //TODO: Figure out the grid allocation for the 3D case
                    M_u = Eigen::TensorXd();
                    M_u.init_tensor(bb_size_x,bb_size_y,bb_size_z);

                    M_v = Eigen::TensorXd();
                    M_v.init_tensor(bb_size_x,bb_size_y,bb_size_z);

                    M_w = Eigen::TensorXd();
                    M_w.init_tensor(bb_size_x,bb_size_y,bb_size_z);
                }



void init_state_2d(const int bb_size_x,const int bb_size_y, 
                const double grid_interval, const int num_particles, Eigen::MatrixXd &M_particles,
                Eigen::MatrixXd &M_p, Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v){
                    M_particles = Eigen::MatrixXd(num_particles,3);

                    M_p = Eigen::MatrixXd(bb_size_x,bb_size_y);
                    M_u = Eigen::MatrixXd(bb_size_x,bb_size_y + 1);
                    M_v = Eigen::MatrixXd(bb_size_x + 1,bb_size_y);

                    M_p.setZero();
                    M_u.setZero();
                    M_v.setZero();
                }