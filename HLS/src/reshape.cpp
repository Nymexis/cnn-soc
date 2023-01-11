#include <reshape.h>

int reshape(int x_size, int y_size, int z_size, double img_in[MAX_SIZE_FM], double img_out[SIZE_OUT]){
    int i = 0;
    for (int x = 0; x < x_size ; x++){
        for (int y = 0 ; y < y_size ; y++){
            for (int z = 0 ; z < z_size ; z++){
                img_out[i] = img_in[y+x*x_size+z*y_size*x_size];
                i++;
            }
        }
    }
    return 1;
}