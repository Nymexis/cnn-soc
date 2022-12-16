#include <imgRead.h>

int readImg(const char* im_filepath, int img_idx, Image* img) {
	FILE* fp = fopen(im_filepath, "rb");
	if (!fp) {
		return -1;
	}

	fseek(fp, img_idx * IMG_SIZE, SEEK_SET);
	
	fread(&img->label, 1, 1, fp);
	
	for (int color = 0; color < 3; ++color) {
	    for (int l = 0; l < IMG_PIX_SIZE; ++l) {
	        for (int c = 0; c < IMG_PIX_SIZE; ++c) {
	            fread(&img->img[l][c][color], 1, 1, fp);
	        }
	    }	
	}

	fclose(fp);
	return 0;

}

int main(int argc, char** argv) {
	if (argc < 2) {
		printf("Syntax : ./imageRead image_idx\n");
		return 1;
	}

	Image img;
	if (readImg("../data/cifar10_data/cifar-10-batches-bin/test_batch.bin", atoi(argv[1]), &img) != 0) {
    		printf("Failed to read image\n");
    		return 1;
	}
	
	return 0;

}
