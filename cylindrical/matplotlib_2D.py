# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 09:37:38 2022

@author: xudonghan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# 数据读取
r = np.loadtxt('coordinate.txt')
E_z = np.loadtxt('result.txt')


# matplot 绘图设置
fig = plt.figure()
img = fig.subplots()
img.set_title('line current')
img.set_ylim(-300, 300)
img.set_ylabel('$E_z$')
img.set_xlabel('$r$')

# 循环绘图
ims = []
for n in range(0,500):
    im = img.plot(r, E_z[n*400,:], '-r')
    ims.append(im)

# 动画绘制
ani = animation.ArtistAnimation(fig, ims , interval=50)
ani.save('line current.gif')
plt.show()

