"""
.. _gif_movie_example:

Create a GIF Movie
~~~~~~~~~~~~~~~~~~
Generate a moving gif from an active plotter.

.. note::
   Use ``lighting=False`` to reduce the size of the color space to avoid
   "jittery" GIFs, especially for the scalar bar.

"""

import numpy as np
import pyvista as pv

# Read data
x, y = np.loadtxt('coordinate.txt')
Data = np.loadtxt('Ez.txt')

# Create and structured surface
X, Y = np.meshgrid(x, y)
E_z = Data[0:201,:]
grid = pv.StructuredGrid(X, Y, E_z)

# Create a plotter object and set the scalars to the Z height
plotter = pv.Plotter(notebook=False, off_screen=True)
plotter.add_mesh(
    grid,
    scalars=E_z.ravel(),
    lighting=True,
    show_edges=False,
    scalar_bar_args={"title": "E_z"},
    clim=[-3, 3],
    cmap = 'jet',
)

# Open a gif
plotter.open_gif("wave.gif")

pts = grid.points.copy()
# Update E_z and write a frame for each updated position
nframe = 401
for n in range(0,400)[:nframe]:
    E_z = Data[n*201:(n+1)*201]
    pts[:, -1] = E_z.ravel()
    plotter.update_coordinates(pts, render=False)
    plotter.update_scalars(E_z.ravel(), render=False)
    # Write a frame. This triggers a render.
    plotter.write_frame()
# Closes and finalizes movie
plotter.close()
