# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 19:47:03 2022

@author: xudonghan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 数据读取
# x, y = np.loadtxt('coordinate.txt')
E_z = np.loadtxt('Ez.txt')
# X, Y = np.meshgrid(x, y)

fig = plt.figure()
img = fig.subplots()
img.set_title('line current')

# 循环绘图
ims = []
for n in range(0,400):
    plt.axis('off')
    im = plt.imshow(E_z[n*201:(n+1)*201,:], cmap = 'rainbow', vmin = -4, vmax = 4)
    ims.append([im])

# 动画绘制
ani = animation.ArtistAnimation(fig, ims , interval = 50)
ani.save('line current.gif')
