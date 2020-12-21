#include <velocity_extrapolation.h>
#include <utilities.h>
void extrapolate_velocity_2d(Eigen::MatrixXd &M_u, Eigen::MatrixXd &M_v, Eigen::MatrixXd &M_fluid)
{
    int x_len = M_u.cols();
    int y_len = M_v.rows();
    int x_len_non_staggered = x_len - 1;
    int y_len_non_staggered = y_len - 1;

    for (int x = 0; x < x_len_non_staggered; x++)
    {
        for (int y = 0; y < y_len_non_staggered; y++)
        {
            int i_idx, j_idx, center_i, center_j;
            get_matrix_index_2d(x, y, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);

            if (M_fluid(i_idx, j_idx) == 1)
            {
                continue;
            }
            int count = 0;
            //updating v
            get_matrix_index_2d(x, y, x_len_non_staggered, y_len, center_i, center_j);
            double v = 0;
            //v_x,y-1
            if (!on_boundary(y - 1, y_len_non_staggered))
            { //bottom
                get_matrix_index_2d(x, y - 1, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (M_fluid(i_idx, j_idx) == 1)
                { //fluid cell
                    get_matrix_index_2d(x, y - 1, x_len_non_staggered, y_len, i_idx, j_idx);
                    v += M_v(i_idx, j_idx);
                    count++;
                }
            }

            //v_x,y+1
            if (!on_boundary(y + 1, y_len_non_staggered))
            { //top
                get_matrix_index_2d(x, y + 1, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (M_fluid(i_idx, j_idx) == 1)
                { //fluid cell
                    get_matrix_index_2d(x, y + 1, x_len_non_staggered, y_len, i_idx, j_idx);
                    v += M_v(i_idx, j_idx);
                    count++;
                }
            }

            //v_x-1,y
            if (!on_boundary(x - 1, x_len_non_staggered))
            { //left
                get_matrix_index_2d(x - 1, y, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (M_fluid(i_idx, j_idx) == 1)
                { //fluid cell
                    get_matrix_index_2d(x - 1, y, x_len_non_staggered, y_len, i_idx, j_idx);
                    v += M_v(i_idx, j_idx);
                    count++;
                }
            }

            //v_x+1,y
            if (!on_boundary(x + 1, x_len_non_staggered))
            { //right
                get_matrix_index_2d(x + 1, y, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (M_fluid(i_idx, j_idx) == 1)
                { //fluid cell
                    get_matrix_index_2d(x + 1, y, x_len_non_staggered, y_len, i_idx, j_idx);
                    v += M_v(i_idx, j_idx);
                    count++;
                }
            }
            if (count > 0)
            {
                M_v(center_i, center_j) = v / count;
            }

            //updating u

            get_matrix_index_2d(x, y, x_len, y_len_non_staggered, center_i, center_j);
            double u = 0;
            count = 0;
            //u_x,y-1
            if (!on_boundary(y - 1, y_len_non_staggered))
            { //bottom
                get_matrix_index_2d(x, y - 1, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (M_fluid(i_idx, j_idx) == 1)
                { //fluid cell
                    get_matrix_index_2d(x, y - 1, x_len, y_len_non_staggered, i_idx, j_idx);
                    u += M_u(i_idx, j_idx);
                    count++;
                }
            }

            //u_x,y+1
            if (!on_boundary(y + 1, y_len_non_staggered))
            { //top
                get_matrix_index_2d(x, y + 1, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (M_fluid(i_idx, j_idx) == 1)
                { //fluid cell
                    get_matrix_index_2d(x, y + 1, x_len, y_len_non_staggered, i_idx, j_idx);
                    u += M_u(i_idx, j_idx);
                    count++;
                }
            }

            //u_x-1,y
            if (!on_boundary(x - 1, x_len_non_staggered))
            { //left
                get_matrix_index_2d(x - 1, y, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (M_fluid(i_idx, j_idx) == 1)
                { //fluid cell
                    get_matrix_index_2d(x - 1, y, x_len, y_len_non_staggered, i_idx, j_idx);
                    u += M_u(i_idx, j_idx);
                    count++;
                }
            }

            //u_x+1,y
            if (!on_boundary(x + 1, x_len_non_staggered))
            { //right
                get_matrix_index_2d(x + 1, y, x_len_non_staggered, y_len_non_staggered, i_idx, j_idx);
                if (M_fluid(i_idx, j_idx) == 1)
                { //fluid cell
                    get_matrix_index_2d(x + 1, y, x_len, y_len_non_staggered, i_idx, j_idx);
                    u += M_u(i_idx, j_idx);
                    count++;
                }
            }
            if (count > 0)
            {
                M_u(center_i, center_j) = u / count;
            }
        }
    }
}
