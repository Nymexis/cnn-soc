#include <perceptron.h>
#include <iostream>

int main(){
    double img_in[SIZE_IN];
    double img_out[SIZE_OUT];
    double kernels[ROM_SIZE];

    if (perceptron(img_in, kernels, img_out)){
        std::cout << "BRAVO" << std::endl;
    }

    return 0;
}