#!/usr/bin/env/ python3

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import animation as animation


def LoadTxt(filename):
    data = np.loadtxt(filename)
    return data

def plot1d(data, images):
    im = plt.plot(data, color='tab:blue')
    return im

if __name__ == '__main__':
    images = []
    fig = plt.figure()
    for i in range(21):
        filename = str(i) + '.txt'
        data = LoadTxt(filename)
        im = plot1d(data, images)
    #     plt.show()
        images.append(im)
    ani = animation.ArtistAnimation(fig, images, interval=50, blit=True)
    ani.save('RJ.mp4', writer='ffmpeg')
