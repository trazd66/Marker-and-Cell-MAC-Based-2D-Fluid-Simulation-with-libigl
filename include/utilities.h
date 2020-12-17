#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>
#include <iostream>
/***
 * 
 * transates the x,y position from the paper into eigen indexes
 * 
 * e.g.
 *get_matrix_index_2d(x_non_staggered,y-1,x_len_non_staggered,y_len,i_idx,j_idx);
 * 
 ***/
void get_matrix_index_2d(const int x ,const int y, 
                                const int len_x,const int len_y, 
                                int &i, int &j); 