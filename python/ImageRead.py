# coding : UTF-8

from matplotlib.pyplot import imread, imshow, show
import numpy as np
import struct


class ImageRead:
    IMG_SIZE     = 3072 + 1
    COMP_SIZE    = 1024
    IMG_PIX_SIZE = 32
    def __init__(
        self,
        im_filepath
    ):
        self.im_filepath = im_filepath

    def compute(self):
        arr = imread(self.im_filepath)
        return arr

    def readImg(self, img_idx):
        img = []

        with open(self.im_filepath, "rb") as f:
            f.read(img_idx * self.IMG_SIZE) # go to the image at img_idx

            img_label = f.read(1)[0]
            
            img = np.zeros((self.IMG_PIX_SIZE,self.IMG_PIX_SIZE,3))

            for color in range(3):
                for l in range(self.IMG_PIX_SIZE):
                    for c in range(self.IMG_PIX_SIZE):
                  
                        img[l][c][color] = f.read(1)[0]

            return img, img_label      


if __name__ == "__main__":
    imgObj = ImageRead("../data/cifar10_data/cifar-10-batches-bin/test_batch.bin")

    img, label = imgObj.readImg(6070)

    imshow(img)
    show()
