# @Author: Juan Trejo <jtrejo13>
# @Date:   2018-04-14T16:24:34-05:00
# @Last modified by:   jtrj13
# @Last modified time: 2018-04-21T15:31:00-05:00

import numpy as np
from matplotlib import pyplot
import pandas as pd
from collections import namedtuple

point = namedtuple('Point', ['x', 'y'])

# ----------
# Grid Setup
# ----------
NDIM = 2

NX = 51
NY = 51
dx = 2 / (NX - 1)
dy = 2 / (NY - 1)
x = np.linspace(0, 2, NX)
y = np.linspace(0, 2, NY)

dt = 0.01

u0 = np.ones((NX, NY))
u1 = np.ones((NX, NY))

s0 = np.ones((NX, NY))
s1 = np.ones((NX, NY))

# -------------------
# Physical Properties
# -------------------

visc  = 1.0
diff  = 1.0
diss  = 1.0

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
        source = np.ones((NX, NY))
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

def dens_step(s1, s0, u, v, source, diff, diss, dt):
    add_source(s0, source, dt)
    diffuse(s0, s1, diff, dt)
    advect(s1, s0, u, v, dt)
    # dissipate(s1, s0, diss, dt)


# ----------
# Add Source
# ----------

def add_source(s, source, dt):
    s = dt*source


# ---------------
# Diffuse Density
# ---------------
def diffuse(s0, s1, diff, dt):
    a = dt*diff*np.size(s0)

    for _ in range(21):
        s1[1:-1, 1:-1] = 1 / (1 + 4 * a) * (s0[1:-1, 1:-1] \
                         + a * (s1[0:-2, 1:-1]+ s1[2:, 1:-1]
                         + s1[1:-1, 0:-2] + s1[1:-1, 0:-2]))
        set_boundary() # TODO: Implement for density


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
    i1 = int(clmap_pos(i0 + 1, 0, NX - 1))
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
