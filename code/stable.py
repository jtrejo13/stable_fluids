# @Author: Juan Trejo <jtrejo13>
# @Date:   2018-04-14T16:24:34-05:00
# @Last modified by:   jtrj13
# @Last modified time: 2018-04-23T12:20:38-05:00

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

NX = 41
NY = 41
dx = 2 / (NX - 1)
dy = 2 / (NY - 1)
x = np.linspace(0, 2, NX)
y = np.linspace(0, 2, NY)
X, Y = np.meshgrid(x, y)

dt = 0.01

u0 = np.zeros((NX, NY))
u1 = np.empty_like(u0)
v0 = np.zeros((NX, NY))
v1 = np.empty_like(u0)

s0 = np.zeros((NX, NY))
s1 = np.empty_like(u0)

# -------------------
# Physical Properties
# -------------------
visc  = 1.0
diff  = 1.0
diss  = 1.0

# ------------------
# Initial Conditions
# ------------------

# Density
s0[int(0.5 / dx):int(0.7 / dx + 1), int(0.5 / dy):int(0.7 / dy + 1)] = 2
s1 = np.copy(s0)
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
        force = np.ones((NX, NY))
        source = np.zeros((NX, NY))
        vel_step(u1, u0, visc, force, dt)
        dens_step(s1, s0, u1, v1, source, diff, diss, dt)

# -------------
# Velocity Step
# -------------
def vel_step(u1, u0, visc, force, dt):
    pass

# ------------
# Density Step
# ------------
def dens_step(s1, s0, u, v, source, diff, diss, dt):
    add_source(s0, source, dt)
    diffuse(s0, s1, diff, dt)
    advect(s1, s0, u, v, dt)
    # dissipate(s1, s0, diss, dt)


# ----------
# Add Source
# ----------
def add_source(s, source, dt):
    s = s + dt * source


# ---------------
# Diffuse Density
# ---------------
def diffuse(s0, s1, diff, dt):
    a = dt * diff * np.size(s0)

    for _ in range(21):
        s1[1:-1, 1:-1] = 1 / (1 + 4 * a) * (s0[1:-1, 1:-1] \
                         + a * (s1[0:-2, 1:-1]+ s1[2:, 1:-1]
                         + s1[1:-1, 0:-2] + s1[1:-1, 0:-2]))
        set_boundary(s1) # TODO: Implement for density


# --------------
# Advect Density
# --------------
def advect(s1, s0, u, v, dt):
    dt0 = -dt * NX

    for i in range(NX):
        for j in range(NY):
            orig_pos = point(i * dx, j * dy)
            orig_vel = point(u[i, j], v[i, j])

            # RK2
            half_dt = 0.5 * dt0
            halfway_pos = point(orig_pos.x + half_dt * orig_vel.x, \
                                orig_pos.y + half_dt * orig_vel.y)

            # Restrict pos values if boundary doesn't wrap around
            halfway_pos.x = clamp_pos(halfway_pos.x, 0, NX - 1)
            halfway_pos.x = clamp_pos(halfway_pos.y, 0, NY - 1)

            halfway_vel = point(interpolate(u, halfway_pos), \
                                interpolate(v, halfway_pos))

            backtrack_pos = point(orig_pos.x + dt0 * halfway_vel.x, \
                                  orig_pos.y + dt0 * halfway_vel.y)

            # Restrict pos values if boundary doesn't wrap around
            backtrack_pos.x = clamp_pos(backtrack_pos.x, 0, NX - 1)
            backtrack_pos.y = clamp_pos(backtrack_pos.y, 0, NY - 1)

            traced_s = interpolate(s0, backtrack_pos)

            s1[i, j] = traced_s

# ----------------
# Helper Functions
# ----------------
def set_boundary(x):
    x[0, :] = 0
    x[:, 0] = 0
    x[:, -1] = 0
    x[-1, :] = 0

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
    return max if pos > max else min if pos < min else pos

def interpolate(arry, pos):
    i0 = int(pos.x / dx)
    j0 = int(pos.y / dy)
    i1 = int(clamp_pos(i0 + 1, 0, NX - 1))
    j1 = int(clamp_pos(j0 + 1, 0, NY - 1))

    it = pos.x - i0 * dx
    jt = pos.y - j0 * dy

    assert (it < 1.0 and it >= 0)
    assert (jt < 1.0 and jt >= 0)

    xBottomLeft   = x[i0, j0]
    xBottomRight  = x[i1, j0]
    xBottomInterp = lerp(xBottomLeft, xBottomRight, it)

    xTopLeft   = x[i0, j1]
    xTopRight  = x[i1, j1]
    xTopInterp = lerp(xTopLeft, xTopRight, it)

    xMidInterp = lerp(xBottomInterp, xTopInterp, jt)

    return xMidInterp

def lerp(v0, v1, t):
    return (1 - t) * v0 + t * v1
