from pprint import pprint
from math import exp

class Perceptron:
    
    input_feature = None
    feature_size  = 0
    biases        = None
    weights       = None
    output_prob   = None
    prob_nb       = 0

    def __init__(
        self,
        input_feature,
        weights,
        biases
        ):
        self.input_feature = input_feature
        self.biases        = biases
        self.weights       = weights
        self.prob_nb       = len(self.biases)
        self.feature_size  = len(self.input_feature)

    def compute(self):
        # print("Perceptron Input")
        # self.print_len(self.input_feature)
        self.output_prob = [0 for i in range(self.prob_nb)]
        for i in range(self.prob_nb):
            for k in range(len(self.input_feature)):
                self.output_prob[i] += self.weights[k][i] * self.input_feature[k]

            #self.output_prob[i] = sum([self.weights[k][i] * self.input_feature[k] \
            #    for k in range(self.prob_nb)])
            self.output_prob[i] += self.biases[i]

        sum_exp = sum([exp(self.output_prob[i]) for i in range(self.prob_nb)])
        for i in range(self.prob_nb):
            self.output_prob[i] = exp(self.output_prob[i]) / sum_exp

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
    
    from Parser import Parser

    # parser = Parser("../data/CNN_coeff_3x3.txt")
    # parser.readCoeff()

    ft = [2.0 for k in range(20)]
    wght = [[1.0 for j in range(20) for k in range(10)]]
    biases = [0. for k in range(10)]

    per = Perceptron(ft, wght, biases) 
    # per = Perceptron(ft, parser.coeffs["local3/weights"], parser.coeffs["local3/biases"])
    per.compute()

    print(per.output_prob)