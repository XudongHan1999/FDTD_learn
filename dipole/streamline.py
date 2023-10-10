# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 11:29:30 2022

@author: xudonghan
"""

# 导入库
import numpy as np
import pyvista as pv

# 电场数据读取
E = np.loadtxt('E.txt')

# 网格节点划分
nx = 60; ny = 60; nz = 60
origin = (-60., -60., -60.)
mesh = pv.UniformGrid(dimensions = (nx, ny, nz), spacing = (2., 2., 2.), origin = origin)

# 电场矢量
vectors = np.empty((mesh.n_points, 3))
for n in range(0, 3):
    vectors[:, n] = E[:, n]
mesh['vectors'] = vectors

# 流场线绘制
stream, src = mesh.streamlines(
    'vectors', return_source=True, terminal_speed=0.0, n_points=200, source_radius=30
)
stream.tube(radius=0.01).plot(cmap = 'jet')
