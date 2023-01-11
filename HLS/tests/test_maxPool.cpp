#include <maxPool.h>
#include <conv2d.h>

int main(){
    int x = 24;
    int y = 24;
    int z = 3;

    double img_in[MAX_F_SIZE];
    double img_out[MAX_F_SIZE];

    if (maxpool(x, y, z, img_in, img_out)){
        std::cout << "BRAVO" << std::endl;
    }  

    return 0;
}