#ifndef _PERCEPTRON_H_
#define _PERCEPTRON_H_

#include <cmath>

const int SIZE_IN  = 180;
const int SIZE_OUT = 10;

const int MAX_SIZE_K = 24*24*64+64;
const int K_NBR_     = 3;
const int END_MATRIX = SIZE_IN*SIZE_OUT+SIZE_OUT;

const int ROM_SIZE = MAX_SIZE_K * K_NBR_ + END_MATRIX;

int perceptron(double img_in[SIZE_IN], double kernels[ROM_SIZE], double img_out[SIZE_OUT]);

#endif