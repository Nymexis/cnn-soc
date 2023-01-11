#ifndef TYPES_H
#define TYPES_H

#include "ac_fixed.h"
#define  WIDTH = 12
#define  I     = 9

// | <------------ W ----------------> |
// | 0 <-- I --> 0 . 0 <-- W - I --> 0 |
typedef ac_fixed<WIDTH, I, true, AC_RND> dType;
// AC_RND -> Round up number towards infinity

// unsigned 64 bits integer
typedef ac_int<6,false> iType;

#endif