# coding : utf-8
'''
Parser class 
read a text file containing coefficients for the convolution operation.
Any coefficient array in the text file must be after a tensor name :
tensor_name:  name
[[v v v v v
  v v v v v]
 [v v v v v
  v v v v v]
    .....
 [v v v v v
  v v v v v]
 [v v v v v
  v v v v v]]

- values in the array must be separated by spaces.
- The array can be written on multiple lines.
- exponential notation is accepted
- new subarray must start on a new line

Parsing result are stored in the coeffs dictionaries with string keys and lists of floats as values
'''
import re
from pprint import pprint

class Parser:

    f_path    = ""
    coeffs = {}

    def __init__(
        self, 
        f_path
        ):
        self.f_path = f_path

    def readCoeff(self):
        f           = open(self.f_path)
        lines       = f.readlines()
        tensor_name = ""
        array_str   = ""

        for l in lines:
            if l.split(':')[0] == "tensor_name":
                if tensor_name != "":
                    self.coeffs[tensor_name] = array_str

                tensor_name = re.sub(" +", '', l.split(":")[1]).replace('\n','')
                array_str = ""
            elif len(l) > 1:
                tmp_line = re.sub("\[ +",'[', l)        # Removing space after/before bracket to avoid wrong comma insertion
                tmp_line = re.sub(" +\]",']', tmp_line) # ...
                tmp_line = re.sub(" +"  ,',', tmp_line).replace('\n', '') # change space to comma
        
                chr_test = re.match("[\[\]0-9\-,.e]+", tmp_line) # sanetization of input

                if chr_test == None or chr_test.string != tmp_line:
                    print(tmp_line)
                    print(chr_test)
                    print("Invalid character detected.") # could also mean empty line
                    return 0

                array_str += tmp_line
        
        self.coeffs[tensor_name] = array_str

        # translate array_str into a real float array
        #istce_name = [k for k,v in globals().items() if v is self][0] # name of the instance variable
        for k,v in self.coeffs.items():
            exec("self.coeffs['"+k+"'] = "+v)
        
        return 1

if __name__ == "__main__":
    from time import time
    from pprint import pprint

    t0 = time()
    prsr = Parser("../data/CNN_coeff_3x3.txt")
    t1 = time()

    print("Total parsing time {} ms".format((t1-t0)*1000))

    if prsr.readCoeff():
        print(prsr.coeffs.keys())
        wght = prsr.coeffs["conv2/weights"]

        pprint(wght)
    else:
        print("An error occured.")
