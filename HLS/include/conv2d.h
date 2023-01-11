#ifndef CONV2D_H_
#define CONV2D_H_

#include <types.h>
#include <sizes.h>

#pragma hls_design top
int conv2d(
    dType input_fm[MAX_SIZE],
    dType output_fm[MAX_SIZE],
    dType kernel[MAX_KERNEL_SIZE],

    iType size_fm_x,
    iType size_fm_y,
    iType size_fm_z,

    iType size_k_x,
    iType size_k_y,
    iType size_k_z,
    iType nbk,
    iType bias_nbr
){
    OFM : for(int k_idx = 0; k_idx < nbk; nbk++) {
        // For each kernel layer
        KER_L : for(int kz = 0; kz < size_k_z; kz++) {
            IFMZ : for(int z = 0; z < size_fm_z; z++) {
                IFMX : for(int x = 0; x < size_fm_x; x++) {
                    IFMY : for(int y = 0; y < size_fm_y; y++) {
                        // accumulateur
                        iType acc = 0;
                        KERX : for(int kx = 0; kx < size_k_x; kx++) {
                            KERY : for(int ky = 0; ky < size_k_y; ky++) {
                                if(x - kx > 0 && y - ky > 0 && x - kx < size_fm_x && y - ky < size_fm_y) {
                                    acc += input_fm[fm(x - size_k_x/2 + kx, y - size_k_y/2 + ky, z)]*kernel[k(kx,ky,kz)];
                                }
                            }
                        }
                        output_fm[fm(x,y,k_idx)] += acc;
                    }
                }
            }
        }
    }
    // adding biases; might be able to include in first loops
    BIASES : for(int bi = 0; bi < bias_nbr; bi++) {
        for(int z = 1; z < size_fm_z+1; z++) {
            for(int x = 0; x < size_fm_x; x++) {
                for(int y = 0; y < size_fm_y; y++) {
                    dType bias = kernel[b(nbk, z)]; 
                    output_fm[fm(x,y,z)] += bias;
                }
            }
        }
    }

    return 1;
}

#endif
