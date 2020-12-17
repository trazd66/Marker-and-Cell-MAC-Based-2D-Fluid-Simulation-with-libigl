#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <EigenTypes.h>

/***
 * 
 * transates the x,y position from the paper into eigen indexes
 * 
 * 
 ***/
void get_matrix_index_2d(const int x ,const int y, 
                                const int len_x,const int len_y, 
                                int &i, int &j); 