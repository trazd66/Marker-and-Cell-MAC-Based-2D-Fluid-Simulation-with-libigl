#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <iostream>

/* adopted from https://stackoverflow.com/questions/3767869/adding-message-to-assert */
#ifndef NDEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

/***
 *
 * transate the x,y position from the paper into eigen indexes
 * e.g. get_matrix_index_2d(x_non_staggered,y-1,x_len_non_staggered,y_len,i_idx,j_idx);
 * Effect: updates i and j.
 ***/
void get_matrix_index_2d(const int x ,const int y, 
                                const int len_x,const int len_y, 
                                int &i, int &j); 


/*
    check if x is out of boundary.
*/
bool on_boundary(const int x, const int len_x);

/*
    check if particle is out of boundary.
*/
bool is_out_of_boundary(Eigen::Vector2d particle_pos, double grid_interval, int len_x, int len_y);

/*
    Ensure particle velocity * dt is strictly less than one grid interval.
    Effect: updates dt.
*/
void update_dt(double &dt, int num_particles, double grid_interval, Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v);

/*
    normalize particle velocities using maximum velocity within each matrix.
    Effect: updates M_particles_u and M_particles_v
*/
void normalize_velocity(Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v);

/*
    normalize grid velocities using maximum velocity within each matrix.
    Effect: updates M_u and M_v.
*/
void normalize_grid(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v);

/*
    helper function to get bilinear interpolation coefficients.
    Effect: updates coeff.
*/
void get_bilinear_coeff(Eigen::Vector4d &coeff, double grid_interval, int grid_x_start, int grid_x_end, int grid_y_start, int grid_y_end, double particle_x, double particle_y);

/*
    get grid index from particle postiion.
    Effect: updates grid_x and grid_y.
*/
void particle_pos_to_grid(int &grid_x, int &grid_y, double x, double y, double grid_interval);