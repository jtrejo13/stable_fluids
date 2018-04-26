# @Author: Juan Trejo <jtrejo13>
# @Date:   2018-04-14T16:24:34-05:00
# @Last modified by:   jtrj13
# @Last modified time: 2018-04-26T13:18:59-05:00

# -------
# Imports
# -------
import numpy as np
import pandas as pd
from collections import namedtuple
from matplotlib import cm, pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

point = namedtuple('Point', ['x', 'y'])

# ----------
# Grid Setup
# ----------
NDIM = 2

NX = 81
NY = 81
dx = 2 / (NX - 1)
dy = 2 / (NY - 1)
x = np.linspace(0, 2, NX)
y = np.linspace(0, 2, NY)
X, Y = np.meshgrid(x, y)

dt = 0.01

# -------------------
# Physical Properties
# -------------------
visc  = 1.0
diff  = 1
diss  = 1.0

# ------------------
# Initial Conditions
# ------------------

u0 = np.zeros((NY, NX))
u1 = np.zeros((NY, NX)) + .1
v0 = np.zeros((NY, NX))
v1 = np.zeros((NY, NX))

s0 = np.zeros((NY, NX))
s1 = np.zeros((NY, NX))
source = np.zeros((NY, NX))

# u1 = np.full((NY, NX), 0.1)
s0[int(0.9 / dy):int(1.1 / dy + 1), int(0.9 / dx):int(1.1 / dx + 1)] = 10
s1[int(0.9 / dy):int(1.1 / dy + 1), int(0.9 / dx):int(1.1 / dx + 1)] = 10
source[int(0.9 / dy):int(1.1 / dy + 1), int(0.1 / dx):int(0.3 / dx + 1)] = 10

# ---------
# Main Loop
# ---------
def solve():
    for i in range(20):
    # Handle display
    # Get forces and sources
    s1, s0 = s0, s1
    # vel_step(u1, u0, visc, force, dt)
    dens_step(s1, s0, u1, v1, source, diff, diss, dt)
    source = np.zeros((NY, NX))

# -------------
# Velocity Step
# -------------
def vel_step(u1, u0, visc, force, dt):
    pass

# ------------
# Density Step
# ------------
def dens_step(s1, s0, u, v, source, diff, diss, dt):
    add_source(s1, s0, dt)
    add_source(s1, source, dt)
    s1, s0 = s0, s1
    diffuse(s1, s0, diff, dt)
    s1, s0 = s0, s1
    advect(s1, s0, u, v, dt)
#     dissipate(s1, s0, diss, dt)

# ----------
# Add Source
# ----------
def add_source(s, source, dt):
    s += dt * source

# ---------------
# Diffuse Density
# ---------------
def diffuse(s1, s0, diff, dt):
    a = dt * diff * np.size(s0)

    for _ in range(20):
        s1[1:-1, 1:-1] = (s0[1:-1, 1:-1] + a * (s1[1:-1, 0:-2] + s1[1:-1, 2:] + s1[0:-2, 1:-1] + s1[2:, 1:-1])) / (1 + 4 * a)
        set_boundary(s1)


# --------------
# Advect Density
# --------------
def advect(s1, s0, u, v, dt):
    dt0 = -dt * NX

    for i in range(1, NX-1):
        for j in range(1, NY-1):
            orig_pos = point(i * dx, j * dy)
            orig_vel = point(u[j, i], v[j, i])

            # RK2
            half_dt = 0.5 * dt0
            halfway_pos_x = orig_pos.x + half_dt * orig_vel.x
            halfway_pos_y = orig_pos.y + half_dt * orig_vel.y

            # Restrict pos values if boundary doesn't wrap around
            halfway_pos = point(clamp_pos(halfway_pos_x, 0, NX - 1), \
                                clamp_pos(halfway_pos_y, 0, NY - 1))

            halfway_vel = point(interpolate(u, halfway_pos), \
                                interpolate(v, halfway_pos))

            backtrack_pos_x = orig_pos.x + dt0 * halfway_vel.x
            backtrack_pos_y = orig_pos.y + dt0 * halfway_vel.y

            # Restrict pos values if boundary doesn't wrap around
            backtrack_pos = point(clamp_pos(backtrack_pos_x, 0, NX - 1), \
                                  clamp_pos(backtrack_pos_y, 0, NY - 1))

            traced_s = interpolate(s0, backtrack_pos)
            s1[j, i] = traced_s

    set_boundary(s1)

# ----------------
# Helper Functions
# ----------------
def set_boundary(x, boundary=''):
    x[0, 1:-1] = -x[1,  1:-1] if boundary == 'y_vel' else x[1,  1:-1]
    x[-1,1:-1] = -x[-2, 1:-1] if boundary == 'y_vel' else x[-2, 1:-1]
    x[1:-1, 0] = -x[1:-1,  1] if boundary == 'x_vel' else x[1:-1,  1]
    x[1:-1,-1] = -x[1:-1, -2] if boundary == 'x_vel' else x[1:-1, -2]

    x[0,   0] = 0.5 * (x[1, 0]   + x[0, 1])
    x[0,  -1] = 0.5 * (x[1, -1]  + x[0, -2])
    x[-1,  0] = 0.5 * (x[-2, 0]  + x[-1, 1])
    x[-1, -1] = 0.5 * (x[-2, -1] + x[-1, -2])

def clamp_pos(pos, minval, maxval):
    """
    Clamp position coordinate acording to min and max limits.

    Parameters
    ----------
    pos : float
        Position coordinate.
    minval : float
        Minimum value allowed.
    maxval : float
        Maximum value allowed.

    Returns
    -------
    float
        Clamped position coordinate.
    """
    return maxval if pos > maxval else minval if pos < minval else pos

def interpolate(x, pos):
    i0 = int(pos.x / dx)
    j0 = int(pos.y / dy)
    i1 = int(clamp_pos(i0 + 1, 0, NX - 1))
    j1 = int(clamp_pos(j0 + 1, 0, NY - 1))

    it = pos.x - i0 * dx
    jt = pos.y - j0 * dy

    assert (it < 1.0 and it >= 0)
    assert (jt < 1.0 and jt >= 0)

    xBottomLeft   = x[j1, i0]
    xBottomRight  = x[j1, i1]
    xBottomInterp = lerp(xBottomLeft, xBottomRight, it)

    xTopLeft   = x[j0, i0]
    xTopRight  = x[j0, i1]
    xTopInterp = lerp(xTopLeft, xTopRight, it)

    xMidInterp = lerp(xBottomInterp, xTopInterp, jt)

    return xMidInterp

def lerp(v0, v1, t):
    return (1 - t) * v0 + t * v1
