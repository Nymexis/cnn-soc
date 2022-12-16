#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define IMG_SIZE 3073
#define COMP_SIZE 1024
#define IMG_PIX_SIZE 32

typedef struct Image {
	uint8_t img[IMG_PIX_SIZE][IMG_PIX_SIZE][3];
	uint8_t label;
} Image;


int readImg(const char* im_filepath, int img_idx, Image* img);