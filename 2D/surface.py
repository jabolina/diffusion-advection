from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

with open('numeric.dat') as f:
    print('[+] Opening file.')
    lines = f.readlines()

    i = 0
    x = []
    y = []
    z = []

    for line in lines:
        xv = np.float128(line.split()[0])
        yv = np.float128(line.split()[1])
        zv = np.float128(line.split()[2])

        x.append(xv)
        y.append(yv)
        z.append(zv)

        if xv == max(y) and xv == yv and xv != 0:
            print('[+] Generating image ' + str(i+1))

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            surf = ax.plot_trisurf(x, y, z, cmap=cm.coolwarm, linewidth=0, vmax=0.9, vmin=0.0)
            plt.savefig('plot/fig' + str(i) + '.png', dpi=270, format='png')
            i += 1

            x = []
            y = []
            z = []

