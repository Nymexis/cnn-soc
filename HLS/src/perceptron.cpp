#include <perceptron.h>
#include <iostream>

int perceptron(double img_in[MAX_SIZE], double kernels[MAX_KERNEL_SIZE], double img_out[MAX_SIZE]){
    int    cursor  = MAX_SIZE_K * K_NBR_;
    double sum_exp = 0;
    for (int i = 0; i < SIZE_OUT ; i++){
        img_out[i] = 0;
        for (int j = 0 ; j < SIZE_IN ; j++){
            img_out[i] += img_in[j] * kernels[MAX_SIZE_K*K_NBR_+i*SIZE_IN+j];
        }
        img_out[i] += kernels[MAX_SIZE_K*K_NBR_+END_MATRIX+i];
        sum_exp += exp(img_out[i]);
    }
    for (int i = 0; i < SIZE_OUT ; i++){
        img_out[i] = exp(img_out[i]) / sum_exp;
    }

    return 1;
}