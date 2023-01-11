#ifndef CONV2D_H_
#define CONV2D_H_

#include "types.h"

#define MAX_SIZE_X 24
#define MAX_SIZE_Y 24
#define MAX_SIZE_Z 64
#define MAX_SIZE MAX_SIZE_X*MAX_SIZE_Y*MAX_SIZE_Z

#define MAX_K_SIZE_X 24
#define MAX_K_SIZE_Y 24
#define MAX_K_SIZE_Z 64
#define MAX_BIAS_NBR 64
#define K_NBR        64
#define END_MATRIX_S 180*10 + 10
#define MAX_KERNEL_SIZE (MAX_K_SIZE_X*MAX_K_SIZE_Y*MAX_K_SIZE_Z+MAX_BIAS_NBR)*K_NBR + END_MATRIX_S

#define fm(x,y,z)   x + y*MAX_K_SIZE_X + z*MAX_K_SIZE_Y*MAX_K_SIZE_X
#define k(x,y,z)    x + y*MAX_K_SIZE_X + z*(MAX_K_SIZE_Y*MAX_K_SIZE_X+MAX_BIAS_NBR)
#define b(k,bi)     (k+1)*(MAX_K_SIZE_X*MAX_K_SIZE_Y*MAX_K_SIZE_Z) + bi  

#pragma hls_design top
int conv2d(
    dType input_fm[MAX_SIZE],
    dType output_fm[MAX_SIZE],
    dType kernel[MAX_KERNEL_SIZE],

    iType size_fm_x,
    iType size_fm_y,
    iType size_fm_z,
    iType output_fm_layer,

    iType size_k_x,
    iType size_k_y,
    iType size_k_z,
    iType nbk,
    iType bias_nbr
){
    OFM : for(int k_idx = 0; k_idx < nbk; nbk++) {
        // For each kernel layer
        KER_L : for(int kl = 0; kl < size_k_z; kl++) {
            IFMZ : for(int z = 0; z < size_fm_z; z++) {
                IFMX : for(int x = 0; x < size_fm_x; x++) {
                    IFMY : for(int y = 0; y < size_fm_y; y++) {
                        // accumulateur
                        int acc = 0;
                        KERX : for(int kx = 0; kx < size_k_x; kx++) {
                            KERY : for(int ky = 0; ky < size_k_y; ky++) {
                                if(x - kx > 0 && y - ky > 0 && x - kx < size_fm_x && y - ky < size_fm_y) {
                                    acc += input_fm[fm(x - size_k_x/2 + kx, y - size_k_y/2 + ky, z)]*kernel[k(kx,ky,kz)];
                                }
                            }
                        }
                        output_fm[fm(x,y,z)] += acc;
                    }
                }
            }
        }
    }
    // adding biases; might be able to include in first loops
    BIASES : for(int bi = 0; bi < bias_nbr; bi++) {
        for(z = 0; z < size_fm_z; z++) {
            for(x = 0; x < size_fm_x, x++) {
                for(y = 0; y < size_fm_y; y++) {
                    output_fm[fm(x,y,z)] += kernel[b(nbk, z)];       
                }
            }
        }
    }
}

#endif
