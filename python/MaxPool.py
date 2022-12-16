# coding: UTF-8

from ImageRead import ImageRead
from pprint import pprint 
import numpy
from matplotlib.pyplot import imread, subplots, show

# Matrice used to have a 3x3 MaxPool
m = [
    [-1, -1],
    [0, -1],
    [1, -1],
    [-1, 0],
    [0, 0],
    [1, 0],
    [-1, 1],
    [0, 1],
    [1, 1]
]

class MaxPool:

    def __init__(
        self,
        input_img,
    ):
        # Input image charac
        self.input_img = input_img
        self.size_x = len(input_img[0])
        self.size_y = len(input_img[0][0])
        self.size_z = len(input_img)

        # Image padding
        self.input_img_ad = [[[0 for i in range(self.size_y+1)] for j in range(self.size_x+1)] for k in range(self.size_z)]
        for z in range(self.size_z):
            for y in range(self.size_y):
                for x in range(self.size_x):
                    self.input_img_ad[z][y][x] = self.input_img[z][y][x]

        #output image charac
        self.output_size_x = int(self.size_x/2)
        self.output_size_y = int(self.size_y/2)
        self.output_volume = [[[float(0) for x in range(int(self.size_x/2))] for y in range(int(self.size_y/2))] for z in range(int(self.size_z))]

    def compute(self):
        z = 0
        while z < self.size_z:
            y = 1
            while y < self.size_y:
                x = 1
                while x < self.size_x:
                    max = 0
                    for e in m:
                        cursor_x = e[1] + x
                        cursor_y = e[0] + y
                        if max < self.input_img_ad[z][cursor_y][cursor_x]:
                            max = self.input_img_ad[z][cursor_y][cursor_x]
                    self.output_volume[z][int(y/2)][int(x/2)] = max
                    x += 2
                y += 2
            z += 1
        return 1
        
if __name__ == "__main__":
    from matplotlib.pyplot import imread, subplots, show

    img = ImageRead("../data/cifar10-voilier.png")
    img_l = img.compute()
    mp = MaxPool(img_l)
    mp.compute()
    mp_l = mp.output_volume
