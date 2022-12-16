#!/usr/bin/python3
# coding : utf-8
# main program for this basic CNN processing

import sys
from pprint import pprint

from Adaptator import Adaptator
from ConvolutionNew import Convolution
from ImageRead import ImageRead
from MaxPool import MaxPool
from Parser import Parser
from Perceptron import Perceptron
from Reshape import Reshape


def script(cnn, img_raw):
    # Adaptating and normalizing image
    adaptator = Adaptator(24, 24)
    adaptator.imgLoad(img_raw)
    adaptator.normalize()
    
    # Convolution #1
    img_conv1 = Convolution(adaptator.output_img, cnn.coeffs["conv1/weights"], cnn.coeffs["conv1/biases"], 0)
    img_conv1.convertKernel()
    img_conv1.extractLayers()
    img_conv1.convolve()
   

    # MaxPool #1
    img_mp1 = MaxPool(img_conv1.output_volume)
    img_mp1.compute()
    
    # Convolution #2
    img_conv2 = Convolution(img_mp1.output_volume, cnn.coeffs["conv2/weights"], cnn.coeffs["conv2/biases"], 0)
    img_conv2.convertKernel()
    img_conv2.convolve()
    
    # MaxPool #2
    img_mp2 = MaxPool(img_conv2.output_volume)
    img_mp2.compute()
    
    # Convolution #3
    img_conv3 = Convolution(img_mp2.output_volume, cnn.coeffs["conv3/weights"], cnn.coeffs["conv3/biases"], 0)
    img_conv3.convertKernel()
    img_conv3.convolve()
    
    # MaxPool #3
    img_mp3 = MaxPool(img_conv3.output_volume)
    img_mp3.compute()
    
    # Reshape
    img_reshape = Reshape(img_mp3.output_volume)
    img_reshape.compute()

    pprint(img_reshape.output_volume)
    # Perceptron
    img_percep = Perceptron(img_reshape.output_volume, cnn.coeffs["local3/weights"], cnn.coeffs["local3/biases"])
    img_percep.output_prob = []
    img_percep.compute()
    
    # Checking and writing result
    class_index = img_percep.output_prob.index(max(img_percep.output_prob))
    
    return class_index
if __name__ == "__main__":
    succeeded = 0
    try:
        iterate = int(sys.argv[1])
    except:
        iterate = 10
    print("# Number of iterations : ", iterate)

    # Parsing of the CNN
    cnn = Parser("../data/CNN_coeff_3x3.txt")
    cnn.readCoeff()
    if not cnn:
        print("An error occured while reading coeffs.")
    
    # Reading of input images
    imgs = ImageRead("../data/cifar-10-batches-bin/data_batch_2.bin")

    # Opening results.txt file to store data
    fp_res = open("./results.txt", "w")

    # Opening and reading batches.meta.txt
    fp_meta = open("../data/cifar10_data/cifar-10-batches-bin/batches.meta.txt", "r")
    fp_meta_l = fp_meta.read().splitlines()

    for index in range(iterate):
        print("#", index+1, "/", iterate)
        img_raw, img_label = imgs.readImg(index)
        class_index = script(cnn, img_raw)
        res = 0
        if class_index == img_label:
            res = 1
            succeeded += 1
        res_str = fp_meta_l[img_label] + ":" + fp_meta_l[class_index] + ":" + str(res) + "\n"
        fp_res.write(res_str)
    print("####  ", 100*(succeeded/iterate), "%  ####")
    
