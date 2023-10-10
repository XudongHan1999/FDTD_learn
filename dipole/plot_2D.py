# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 10:24:04 2022

@author: xudonghan
"""

# 导入库
import numpy as np
import matplotlib.pyplot as plt

# 数据读取
x, y, z = np.loadtxt('coordinate.txt')
Y, Z = np.meshgrid(y, z)
E_y = np.loadtxt('E_y.txt')
E_z = np.loadtxt('E_z.txt')

# matplot 绘图设置
fig = plt.figure()
img = fig.subplots()
img.set_title('Hertz Diople')
img.set_xlabel('$y$'); img.set_ylabel('$z$')
img.set_aspect(1)

# 绘制流场线
img.streamplot(Y, Z, E_y, E_z, density = 2.5, linewidth = 0.7, color = 'r')
plt.show()

fig.savefig('dipole.svg', transparent=True)