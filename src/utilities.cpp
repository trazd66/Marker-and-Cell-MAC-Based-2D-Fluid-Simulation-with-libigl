#include <utilities.h>

void get_matrix_index_2d(const int x ,const int y, 
                                const int len_x,const int len_y, 
                                int &i, int &j){
                                    i = len_y - y - 1;
                                    j = x;
                                }


bool on_boundary(const int x, const int len_x){
    return (x >= len_x-1 || x <= 0);
}