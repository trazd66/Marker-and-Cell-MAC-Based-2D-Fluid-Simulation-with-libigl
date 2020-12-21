#include <utilities.h>
#include <algorithm>
#include <iostream>

void get_matrix_index_2d(const int x, const int y,
						 const int len_x, const int len_y,
						 int &i, int &j)
{
	i = len_y - y - 1;
	j = x;
}

bool on_boundary(const int x, const int len_x)
{
	return (x >= len_x - 1 || x <= 0);
}

bool is_out_of_boundary(Eigen::Vector2d particle_pos, double grid_interval, int len_x, int len_y)
{
	double particle_x = particle_pos[0];
	double particle_y = particle_pos[1];
	return particle_x < 0 || particle_x > grid_interval * len_x || particle_y < 0 || particle_y > grid_interval * len_y;
}

/* Ensures particle velocity * dt is strictly less than one grid interval */
void update_dt(double &dt, int num_particles, double grid_interval, Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v)
{
	double u_max = M_particles_u.cwiseAbs().maxCoeff();
	double v_max = M_particles_v.cwiseAbs().maxCoeff();

	double dt_u_tmp;
	if (u_max == 0)
	{
		dt_u_tmp = dt;
	}
	else
	{
		dt_u_tmp = grid_interval / u_max;
		if (dt_u_tmp > dt)
		{
			dt_u_tmp = dt;
		}
	}

	double dt_v_tmp;
	if (v_max == 0)
	{
		dt_v_tmp = dt;
	}
	else
	{
		dt_v_tmp = grid_interval / v_max;
		if (dt_v_tmp > dt)
		{
			dt_v_tmp = dt;
		}
	}

	dt = std::min(dt_u_tmp, dt_v_tmp);
}

void normalize_velocity(Eigen::VectorXd &M_particles_u, Eigen::VectorXd &M_particles_v)
{
	double u_max = M_particles_u.cwiseAbs().maxCoeff();
	double v_max = M_particles_v.cwiseAbs().maxCoeff();

	if (u_max != 0)
	{
		M_particles_u /= u_max;
	}

	if (v_max != 0)
	{
		M_particles_v /= v_max;
	}
}

void normalize_grid(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v)
{
	double u_max = M_u.cwiseAbs().maxCoeff();
	double v_max = M_v.cwiseAbs().maxCoeff();

	if (u_max != 0)
	{
		M_u /= u_max;
	}

	if (v_max != 0)
	{
		M_v /= v_max;
	}
}

void get_bilinear_coeff(Eigen::Vector4d &coeff, double grid_interval, int grid_x_start, int grid_x_end, int grid_y_start, int grid_y_end, double particle_x, double particle_y)
{
	// for the sake of readability
	double x_start_y_start_portion = (1 / (grid_interval * grid_interval)) * (grid_x_end * grid_interval - particle_x) * (grid_y_end * grid_interval - particle_y);
	double x_start_y_end_portion = (1 / (grid_interval * grid_interval)) * (grid_x_end * grid_interval - particle_x) * (particle_y - grid_y_start * grid_interval);
	double x_end_y_start_portion = (1 / (grid_interval * grid_interval)) * (particle_x - grid_x_start * grid_interval) * (grid_y_end * grid_interval - particle_y);
	double x_end_y_end_portion = (1 / (grid_interval * grid_interval)) * (particle_x - grid_x_start * grid_interval) * (particle_y - grid_y_start * grid_interval);

	coeff[0] = x_start_y_start_portion;
	coeff[1] = x_start_y_end_portion;
	coeff[2] = x_end_y_start_portion;
	coeff[3] = x_end_y_end_portion;
}

void particle_pos_to_grid(int &grid_x, int &grid_y, double x, double y, double grid_interval){
	grid_x = (int)(x / grid_interval);
	grid_y = (int)(y / grid_interval);

}
