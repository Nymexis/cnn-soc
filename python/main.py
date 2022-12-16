#!/usr/bin/python3
# coding : utf-8
# main program for this basic CNN processing

from ImageRead import ImageRead
from Adaptator import Adaptator
from ConvolutionNew import Convolution
from Parser import Parser
from MaxPool import MaxPool
from Reshape import Reshape
from Perceptron import Perceptron
from matplotlib.pyplot import imread, subplots, show, imshow
import matplotlib.pyplot as plt
from ImageRead import ImageRead
from random import randint
import sys

debug = 1

def main():
    try:
        index = int(sys.argv[1])
    except:
        # index = 0
        index = randint(0,9999)

    from pprint import pprint
    prsr = Parser("../data/CNN_coeff_3x3.txt")

    if not prsr.readCoeff():
        print("An error occured while reading coeffs.")
    
    imgs = ImageRead("../data/cifar-10-batches-bin/data_batch_1.bin")
    img_raw, img_label = imgs.readImg(index)

    # print(img_label)

    adaptator = Adaptator(24, 24)
    adaptator.imgLoad(img_raw)
    adaptator.normalize()

    # print("***** CONV 1 *****\n")
    img_conv1 = Convolution(adaptator.output_img, prsr.coeffs["conv1/weights"], prsr.coeffs["conv1/biases"], 0)
    # print("Number of output layers : "+str(len(prsr.coeffs["conv1/weights"][0][0][0])))
    img_conv1.convertKernel()
    # print("Converted Kernel len    : Canaux : {} ; sX : {} ; sY : {} ".format(len(img_conv1.kernels[0]), len(img_conv1.kernels[0][0]), len(img_conv1.kernels[0][0][0])))
    img_conv1.extractLayers()
    img_conv1.convolve()
        
    img_mp1 = MaxPool(img_conv1.output_volume)
    img_mp1.compute()
    fig, ax = plt.subplots(1, 10, figsize=(7,4))
    for i in range(10):
        im = ax[i].imshow(img_mp1.output_volume[i])
    plt.show()
    # print("***** CONV 2 *****\n")
    img_conv2 = Convolution(img_mp1.output_volume, prsr.coeffs["conv2/weights"], prsr.coeffs["conv2/biases"], 0)
    # print("Number of output layers : "+str(len(prsr.coeffs["conv2/weights"][0][0][0])))
    img_conv2.convertKernel()
    # print("Converted Kernel len    : Canaux : {} ; sX : {} ; sY : {} ".format(len(img_conv2.kernels[0]), len(img_conv2.kernels[0][0]), len(img_conv2.kernels[0][0][0])))
    img_conv2.convolve()
    # print("Output res sX : {} ; sY : {} ; sZ : {}".format(len(img_conv2.output_volume), len(img_conv2.output_volume[0]), len(img_conv2.output_volume[0][0])))
    img_mp2 = MaxPool(img_conv2.output_volume)
    img_mp2.compute()

    # print("***** CONV 3 *****\n")
    img_conv3 = Convolution(img_mp2.output_volume, prsr.coeffs["conv3/weights"], prsr.coeffs["conv3/biases"], 0)
    img_conv3.convertKernel()
    img_conv3.convolve()
    
    img_mp3 = MaxPool(img_conv3.output_volume)
    img_mp3.compute()

    img_reshape = Reshape(img_mp3.output_volume)
    img_reshape.compute()

    img_percep = Perceptron(img_reshape.output_volume, prsr.coeffs["local3/weights"], prsr.coeffs["local3/biases"])
    img_percep.compute()




    # print("****** RESULT ****** \n\n")
    # pprint(img_percep.output_prob)
    index = img_percep.output_prob.index(max(img_percep.output_prob))
    # print("Max is : ", img_percep.output_prob.index(max(img_percep.output_prob)))
    fp = open("../data/cifar10_data/cifar-10-batches-bin/batches.meta.txt")
    fp_l = fp.read().splitlines()

    # print(fp_l[img_label], "->  ", fp_l[index])

    # fig, axs = subplots(1, 2, figsize=(7, 4))
    # fig.canvas.manager.set_window_title('Convolution result')
    fp = open("./results.txt", "a")
    res = str(int(index == img_label))+"\n"
    fp.write(res+ " Found " + fp_l[index] + " was : " + fp_l[img_label])
    if debug == 1:
        plt.imshow(img_raw/255.)
        plt.title(fp_l[index])
        plt.show()
    
    # pprint(img_percep.output_prob)

def print_len(l):
    print("d1 = ", len(l))
    print("d2 = ", len(l[0]))
    print("d3 = ", len(l[0][0]))
    try:
        print("d4 = ", len(l[0][0][0]))
    except:
        print("d4 = None")

if __name__ == "__main__":
    if debug == 1:
        while 1:
            main()
    else:
        main()
