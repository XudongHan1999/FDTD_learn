import numpy as np
import matplotlib.pyplot as plt

Ez = np.loadtxt('Ez.txt')

fig = plt.figure()
img = fig.subplots()
# img.set_title('metal_scatter')

img.imshow(Ez, cmap = 'jet')

plt.axis('off')


plt.show()
fig.savefig('scatter.svg', transparent = 'True')