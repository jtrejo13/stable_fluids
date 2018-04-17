# @Author: Juan Trejo <jtrejo13>
# @Date:   2018-04-14T16:24:34-05:00
# @Last modified by:   jtrejo13
# @Last modified time: 2018-04-16T12:36:49-05:00

import numpy as np
from matplotlib import pyplot
import pandas as pd

# ----------
# Grid Setup
# ----------
ndim = 2

nx = 51
ny = 51
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

dt = 0.01

u0 = np.ones((nx, ny))
u1 = np.ones((nx, ny))

s0 = np.ones((nx, ny))
s1 = np.ones((nx, ny))

# -------------------
# Physical Properties
# -------------------

visc  =
diff  =
diss  =

# -------------------
# Boundary Conditions
# -------------------




# ---------
# Main Loop
# ---------

def solve():
    while(True):
        # Handle display
        # Get forces and sources
        # Swap U and S
        force = np.ones((nx, ny))
        source = np.ones((nx, ny))
        vel_step(u1, u0, visc, force, dt)
        dens_step(s1, s0, diff, diss, u1, source, dt)

# -------------
# Velocity Step
# -------------

def vel_step(u1, u0, visc, force, dt):
    pass

# ------------
# Density Step
# ------------

def dens_step(s1, s0, diff, diss, u, source, dt):
    add_source(s0, source, dt)
    transport(s1, s0, u, dt)
    diffuse(s0, s1, diff, dt)
    # dissipate(s1, s0, diss, dt)


# ---------
# Add Source
# ---------

def add_source(s, source, dt):
    s = dt*source


# -------
# Diffuse
# -------

def diffuse(s0, s1, diff, dt):
    a = dt*diff*np.size(s0)

    for _ in range(21):
        s1[1:-1, 1:-1] = 1 / (1 + 4 * a) * (s0[1:-1, 1:-1] \
                         + a * (s1[0:-2, 1:-1]+ s1[2:, 1:-1]
                         + s1[1:-1, 0:-2] + s1[1:-1, 0:-2]))
        set_boundary() # TODO: Implement for density
