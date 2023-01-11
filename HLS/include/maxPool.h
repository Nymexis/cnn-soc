#ifndef _MAXPOOL_H_
#define _MAXPOOL_H_

const int SIZE_X = 24;
const int SIZE_Y = 24;
const int SIZE_Z = 64;

const int MAX_F_SIZE = SIZE_X * SIZE_Y * SIZE_Z;
#include <iostream>

int maxpool(int x, int y, int z, double img_in[MAX_F_SIZE], double img_out[MAX_F_SIZE]);

#endif