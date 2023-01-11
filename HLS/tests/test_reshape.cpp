#include <reshape.h>
#include <iostream>

int main(){
    double img_in[MAX_SIZE_FM];
    double img_out[SIZE_OUT];

    if (reshape(SIZE_X, SIZE_Y, SIZE_Z, img_in, img_out)){
        std::cout << "BRAVO" << std::endl;
    }
}