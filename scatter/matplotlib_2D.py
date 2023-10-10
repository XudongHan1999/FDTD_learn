# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 10:46:25 2022

@author: xudonghan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 数据读取
"""
# x, y = np.loadtxt('coordinate.txt')
E_z = np.loadtxt('Ez.txt')
X, Y = np.meshgrid(x, y)
"""
Ez = np.loadtxt('Ez.txt')

fig = plt.figure()
img = fig.subplots()
img.set_title('medium scatter')

# 循环绘图
ims = []
for n in range(0,799):
    plt.axis('off')
    im = plt.imshow(Ez[n*161:(n+1)*161,:], cmap = 'rainbow', vmin = -300, vmax = 300)
    ims.append([im])
    


# 动画绘制
ani = animation.ArtistAnimation(fig, ims, interval = 25)
ani.save('medium_scatter_2.mp4')
# plt.show()

