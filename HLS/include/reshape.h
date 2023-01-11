#ifndef _RESHAPE_H_
#define _RESHAPE_H_

#include <types.h>
#include <sizes.h>

int reshape(int x, int y, int z, double img_in[MAX_SIZE], double img_out[RESHAPE_OUT]);
 
#endif