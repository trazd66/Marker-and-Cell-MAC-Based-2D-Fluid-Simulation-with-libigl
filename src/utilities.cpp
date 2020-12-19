#include <utilities.h>
#include <algorithm>
#include <iostream>

void get_matrix_index_2d(const int x ,const int y, 
                                const int len_x,const int len_y, 
                                int &i, int &j){
                                    i = len_y - y - 1;
                                    j = x;
                                }


bool on_boundary(const int x, const int len_x){
    return (x >= len_x-1 || x <= 0);
}

/* Ensures particle velocity * dt is strictly less than one grid interval */
void update_dt(double &dt, int num_particles, double grid_interval, Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v){
	double u_max = M_particles_u.cwiseAbs().maxCoeff();
	double v_max = M_particles_v.cwiseAbs().maxCoeff();

    std::cout << "u_max -> " << u_max << std::endl;
    std::cout << "v_max -> " << v_max << std::endl;

	double dt_u_tmp;
	if (u_max == 0){
		dt_u_tmp = dt;
	} else {
		dt_u_tmp = grid_interval / u_max;
		if (dt_u_tmp > dt){
			dt_u_tmp = dt;
		}
	}

    double dt_v_tmp;
	if (v_max == 0 ){
		dt_v_tmp = dt;
	} else {
		dt_v_tmp  = grid_interval / v_max;
		if (dt_v_tmp > dt){
			dt_v_tmp = dt;
		}
	}

	dt = std::min(dt_u_tmp, dt_v_tmp);
}

void normalize_velocity(Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v) {
	double u_max = M_particles_u.cwiseAbs().maxCoeff();
	double v_max = M_particles_v.cwiseAbs().maxCoeff();

    if (u_max != 0) {
        M_particles_u /= u_max;
    }

    if (v_max != 0) {
        M_particles_v /= v_max;
    }
}

void normalize_grid(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v) {
	double u_max = M_u.cwiseAbs().maxCoeff();
	double v_max = M_v.cwiseAbs().maxCoeff();

    if (u_max != 0) {
        M_u /= u_max;
    }

    if (v_max != 0) {
        M_v /= v_max;
    }
}
