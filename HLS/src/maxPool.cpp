#include <maxPool.h>
#include <iostream>

int maxpool(int x, int y, int z, double img_in[MAX_F_SIZE], double img_out[MAX_F_SIZE]){
    int cursor_x, cursor_y, max;
    if ( x == 0 || y == 0 || z == 0 || img_in == NULL || img_out == NULL){
        return 0;
    }
    else {  
        for (int f = 0 ; f < z ; f++){
            for (int c = 1 ; c < x ; c+=2){
                for (int l = 1 ; l < y ; l+=2){
                    max = 0;
                    for (int d1 = -1 ; d1 <= 1 ; d1++){
                        for (int d2 = -1 ; d2 <= 1 ; d2++){
                            cursor_x = c + d1;
                            cursor_y = l + d2;
                            if (cursor_x > 0 && cursor_y > 0 && cursor_x < x && cursor_y < y ){
                                std::cout << cursor_x << " " << cursor_y << std::endl;
                                if (img_in[f*SIZE_X*SIZE_Y+cursor_x+cursor_y*SIZE_X] > max){
                                    max = img_in[f*SIZE_X*SIZE_Y+cursor_x+cursor_y*SIZE_X];
                                }
                            }
                        }
                    }
                    img_out[f*SIZE_X*SIZE_Y+cursor_x/2+cursor_y*SIZE_X/2] = max;
                }
            }
        }
        return 1;
    }
}