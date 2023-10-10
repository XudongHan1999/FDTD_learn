# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 22:07:01 2022

@author: xudonghan
"""
# 导入库
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 数据读取
Z = np.loadtxt('coordinate.txt')
Hy = np.loadtxt('Hy.txt')

# matplot 绘图设置
fig = plt.figure()
img = fig.subplots()
img.set_title('incident wave')
img.set_ylim(-1.2, 1.2)
img.set_ylabel('$H_y$')
img.set_xlabel('$z$')

# 循环绘图
ims = []
for n in range(0, 999):
    im = img.plot(Z, Hy[n,:], '-r')
    ims.append(im)

# 动画绘制
ani = animation.ArtistAnimation(fig, ims, interval = 25)
ani.save('incident_1D.mp4')
# plt.show()
