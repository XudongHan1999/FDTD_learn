# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 09:37:38 2022

@author: xudonghan
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# 数据读取
data = np.loadtxt('result.txt')

fig = plt.figure()

# matplot 绘图设置
img = fig.subplots()
img.set_title('surface current')
img.set_ylim(-1.5,1.5)
img.set_ylabel('$H_y$')
img.set_xlabel('$z$')


ims = []
# 循环绘图
for n in range(1,402):
    im = img.plot(data[0,:], data[n,:], '-r')
    ims.append(im)

# 动画绘制
ani = animation.ArtistAnimation(fig, ims , interval=50)
ani.save('surface current.gif')
plt.show()
