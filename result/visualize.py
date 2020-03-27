#!/usr/bin/env/ python3

import numpy as np

from matplotlib import pyplot as plt

def LoadTxt(filename):
    data = np.loadtxt(filename)
    return data

def plot1d(data):
    plt.plot(data)
    plt.show()

if __name__ == '__main__':
    filename = '20.txt'
    data = LoadTxt(filename)
    plot1d(data)
    