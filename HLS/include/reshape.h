#ifndef _RESHAPE_H_
#define _RESHAPE_H_

const int MAX_SIZE_FM = 24*24*64;

const int SIZE_X = 3;
const int SIZE_Y = 3;
const int SIZE_Z = 20;

const int SIZE_OUT = 180;

int reshape(int x, int y, int z, double img_in[MAX_SIZE_FM], double img_out[SIZE_OUT]);
 
#endif