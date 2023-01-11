# coding : utf-8
'''
Note :
- correct padding maintain output and input resolution equals depending on the kernel size and striding

    (I + 2*P - K)/S + 1 = OR 

=> ((OR - 1)*S + K - I)/2 = P

OR : Output Resolution should be an integer
InputImg(24x24) and Kernel(3x3) with Stride(1x1) => P = 1 for OutputImg(24x24)

- convolution
Img*tranpose(K) + bias
'''
import numpy as np
from scipy.signal import convolve2d 
from pprint import pprint
from Adaptator import Adaptator
from MaxPool import MaxPool

class Convolution:

    input_volume  = None
    output_volume = []
    kernels       = None
    bias          = None
    padding       = 0

    def __init__(
            self, 
            input_volume, 
            kernels, 
            bias, 
            padding  # or ask for output resolution and calculate padding corresponding ? 
        ):
        self.input_volume  = input_volume
        self.kernels       = kernels
        self.bias          = bias
        self.padding       = padding

    # checked
    def convertKernel(self):
        size_y    = len(self.kernels)
        size_x    = len(self.kernels[0])
        nb_pixel_c= len(self.kernels[0][0]) # (= RGB)
        nb_kernel = len(self.kernels[0][0][0])
        newKernels= np.zeros((nb_kernel, nb_pixel_c, size_x, size_y))

        for x in range(size_x):
            for y in range(size_y):
                for nbc in range(nb_pixel_c):
                    for kl in range(nb_kernel):
                        newKernels[kl][nbc][x][y] = self.kernels[y][x][nbc][kl]
        
        self.kernels = newKernels

    # checked
    def extractLayers(self):
        # extract every layers
        l_nb   = len(self.input_volume[0][0]) # layers numbers 
        layers = [[] for i in range(l_nb)]
        
        for i in range(l_nb):
            for l in self.input_volume:
                layers[i].append([])
                for c in l:
                    layers[i][-1].append(c[i])
        
        self.input_volume = layers

    def regroupLayers(self):
        # regroup all layers
        x_len = len(self.output_volume[0][0])
        y_len = len(self.output_volume[0])
        z_len = len(self.output_volume)
        
        img = []
        for y in range(y_len):
            img.append([])    
            for x in range(x_len):
                img[-1].append([self.output_volume[i][y][x] for i in range(z_len)])

        self.output_volume = img
        
    # checked
    def convolve(self):
        # print("Convolve Input")
        # self.print_len(self.input_volume)
        bias_idx  = 0
        volume    = []
        for mkernel in self.kernels:
            c_layers = np.zeros((len(self.input_volume[0]),len(self.input_volume[0][0])))
            
            l_idx    = 0
            
            #print(len(mkernel))
                for kernel_layer in mkernel:
                # - Convolution of every volume layers
                tK     = np.transpose(kernel_layer)
                # tK = kernel_layer
                # convolve
                c_layers += convolve2d(self.input_volume[l_idx], tK, boundary="fill", mode="same", fillvalue=0) # to replace with custom
                l_idx += 1

            c_layers += self.bias[bias_idx]
            bias_idx += 1

            # relu
            for i in range(len(c_layers)):
                for j in range(len(c_layers[0])):
                    c_layers[i][j] = c_layers[i][j] if c_layers[i][j] > 0 else 0
        
            volume.append(c_layers)
        self.output_volume = volume

    def print_len(self, l):
        print("d1 = ", len(l))
        print("d2 = ", len(l[0]))
        print("d3 = ", len(l[0][0]))
        try:
            print("d4 = ", len(l[0][0][0]))
        except:
            print("d4 = None")

if __name__ == "__main__":
    from matplotlib.pyplot import imread, subplots, show

    # img_path = "tests/airplane4.png"
    # img      = imread(img_path)
    # ad       = Adaptator(32,32)
    # ad.imgLoad(img)
    # ad.normalize()
    img_test = [[[0., 0., 0.],[0.5, 1.0, 1.0],[0., 0., 0.],[0., 0., 0.]],
          [[0., 0., 0.],[1.0, 1.0, 1.0],[0., 0., 0.],[0., 0., 0.]],
          [[0., 0., 0.],[1.0, 1.0, 1.0],[0., 0., 0.],[0., 0., 0.]],
          [[0., 0., 0.],[1.0, 1.0, 1.0],[0., 0., 0.],[0., 0., 0.]]]
    #img_test = ad.output_img

    k1       = [[ 1./9, 1./9, 1./9], # Blur
                [ 1./9, 1./9, 1./9],
                [ 1./9, 1./9, 1./9]]

    k2       = [[ 0,-1, 0], # Sharpens
                [-1, 5,-1],
                [ 0,-1, 0]] 


    k3       = [[ 0, 0, 0], # Identity
                [ 0, 1, 0],
                [ 0, 0, 0]] 

    
    k4       = [[ 0, 0, 0], # Null
                [ 0, 0, 0],
                [ 0, 0, 0]] 

    k1 = k3

    kernel   = [k3, k3, k4]

    bias     = np.zeros(len(img_test) * len(img_test[0]))
    padding  = 0

    conv = Convolution(img_test, [kernel], bias, padding)
    conv.extractLayers()
    conv.convolve()
    conv.regroupLayers()

    fig, axs = subplots(1, 2, figsize=(7, 4))
    fig.canvas.manager.set_window_title('Convolution result')
    axs[1].imshow(conv.output_volume)
    axs[1].set(title="Output image")
    axs[0].imshow(img_test)
    axs[0].set(title="Input image")
    show()