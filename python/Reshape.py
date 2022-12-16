# coding: UTF-8

from pprint import pprint

class Reshape:

    def __init__(
        self,
        input_img,
    ):
        self.input_img = input_img
        self.output_volume = []

    def compute(self):
        # print("Reshape Input")
        # self.print_len(self.input_img)
        for x in range(len(self.input_img[0][0])):
            for y in range(len(self.input_img[0])):
                for z in range(len(self.input_img)):
                    self.output_volume.append(self.input_img[z][y][x])
        return 1
    
    def print_len(self, l):
        try:
            print("d1 = ", len(l))
        except:
            print("d1 = None")
        try:
            print("d2 = ", len(l[0]))
        except:
            print("d2 = None")
        try:
            print("d3 = ", len(l[0][0]))
        except:
            print("d3 = None")
        try:
            print("d4 = ", len(l[0][0][0]))
        except:
            print("d4 = None")

if __name__ == "__main__":
    img = [[[1., 0., 0.],[2., 1.0, 1.0],[4., 0., 0.]],
          [[5., 0., 0.],[1.0, 1.0, 1.0],[0., 0., 0.]],
          [[6., 0., 0.],[1.0, 1.0, 1.0],[0., 0., 0.]],
          [[7., 0., 0.],[1.0, 1.0, 1.0],[0., 0., 0.]]]

    print("**************************************\n")
    print(img)

    print("**************************************\n")
    reshape = Reshape(img)
    reshape.compute()
    pprint(reshape.output_volume)
