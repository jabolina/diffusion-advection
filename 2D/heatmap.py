import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


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
            xi = np.linspace(min(x), max(x), 1000)
            yi = np.linspace(min(y), max(y), 1000)
            zi = griddata((x, y), z, (xi[None,:], yi[:, None]), method='cubic')

            cs = plt.contourf(xi, yi, zi, 15, cmap=plt.cm.rainbow, vmax=0.9, vmin=0.0)
            #plt.colorbar()
            plt.savefig('images/fig' + str(i) + '.png', dpi=270, format='png')
            i += 1

            x = []
            y = []
            z = []

