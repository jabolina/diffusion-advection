import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


with open('numeric.dat') as f:
    print('[+] Opening file.')
    lines = f.readlines()

    x = []
    y = []
    i = 0
    for line in lines:
        xv = np.float(line.split()[0])
        yv = np.float(line.split()[1])

        if xv in x:
            print('[+] Generating image ' + str(i + 1))
            x = np.array(x)
            y = np.array(y)

            fig, (ax,ax2) = plt.subplots(nrows=2, sharex=True)

            extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
            ax.imshow(y[np.newaxis,:], cmap="plasma", aspect="auto", extent=extent)
            ax.set_yticks([])
            ax.set_xlim(extent[0], extent[1])

            ax2.plot(x,y)

            plt.tight_layout()
            plt.savefig('images/fig' + str(i) + '.png', dpi=270, format='png')
            i += 1

            x = []
            y = []

        x.append(xv)
        y.append(yv)

