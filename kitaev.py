#!/usr/bin/env python

"""Exact eigen-energy of the Kitaev model."""

__author__ = "Satoshi Morita"
__version__ = "2.0.0"

import sys
import numpy as np
import scipy.linalg as sl

L1 = 6 # system length along e1
L2 = 2 # system length along e2
M = 2 # skewness of periodic boundary

Jx = 1.0
Jy = 1.0
Jz = 1.0

class Bond():
    def __init__(self,loop=0):
        """Initialize bond configuration."""
        self.bond = np.array([[Jx,Jy,Jz] for i in range(L1*L2)])
        self.flip_loop(loop)


    def flip_loop(self,loop):
        """Insert Wilson loops."""
        if loop==1:
            self.flip_loop_1()
        elif loop==2:
            self.flip_loop_2()
        elif loop==3:
            self.flip_loop_1()
            self.flip_loop_2()

    def flip_loop_1(self):
        """Insert a Wilson loop along L1*e1."""
        # flip y-bond
        for x in range(L1):
            self.bond[index(x,0)][1] *= -1


    def flip_loop_2(self):
        """Insert a Wilson loop along L2*e2+M*e1."""
        # flip x-bond
        for y in range(L2):
            self.bond[index(0,y)][0] *= -1
        # flip y-bond
        for x in range(M):
            self.bond[index(x+1,L2-1)][1] *= -1


    def create_vortex(self,i,j):
        """Create vortices at hexagones i and j"""
        xi,yi = position(i)
        xj,yj = position(j)
        self.create_vortex_1(xi,xj,yi)
        self.create_vortex_2(xj,yi,yj)


    def create_vortex_1(self,x1,x2,y):
        """Create vortices at (x1,y) and (x2,y)"""
        if x1>x2: x1,x2 = x2,x1
        for x in range(x1,x2):
            self.bond[index(x+1,y)][1] *= -1


    def create_vortex_2(self,x,y1,y2):
        """Create vortices at (x,y1 and (x,y2)"""
        if y1>y2: y1,y2 = y2,y1
        for y in range(y1,y2):
            self.bond[index(x,y+1)][0] *= -1


def index(x, y, dr=(0,0)):
    """Return the index of a unit cell."""
    x0 = x+dr[0]
    y0 = y+dr[1]
    x1 = (x0-M*(y0/L2))%L1
    y1 = y0%L2
    return x1+L1*y1


def position(i, dr=(0,0)):
    """Return the coordinate of a unit cell."""
    x0 = (i%L1)+dr[0]
    y0 = (i/L1)+dr[1]
    return ((x0-M*(y0/L2))%L1, y0%L2)


def nn_1(i):
    """Return the index of the nearest neighbor unit cell along e1."""
    x, y = position(i)
    return index(x+1,y)


def nn_2(i):
    """Return the index of the nearest neighbor unit cell along e2."""
    x, y = position(i)
    return index(x,y+1)


def min_energy(bond):
    """Calculate minimum energy.

    Args:
        bond: an instance of Bond or array[L1*L2][3].
    """
    N_unit = L1*L2
    coupling = bond.bond if isinstance(bond, Bond) else bond

    # Create matrix A
    a = np.zeros((N_unit, N_unit), dtype=float)
    for i in range(N_unit):
        a[i][nn_1(i)] += coupling[i][0]
        a[i][nn_2(i)] += coupling[i][1]
        a[i][i] += coupling[i][2]

    u,s,vt = sl.svd(a)
    det_u = sl.det(u)
    det_vt = sl.det(vt)

    # calculate parity of the projection operator
    ## product of u_{ij}
    sgn = np.prod(np.sign(coupling))
    ## from boundary condition
    if (L1+L2+M*(L1-M))%2 != 0: sgn *= -1 # (-1)^theta
    ## det(Q) = det(VU)
    sgn *= det_u*det_vt

    min_epsilon = min(s)
    sum_epsilon = -0.5*sum(s)
    ene_phys = sum_epsilon
    ene_unphys = sum_epsilon + min_epsilon

    # judge whether the vacuume state is physical or not
    if sgn < 0: # The vacuum state is unphysical.
        ene_phys, ene_unphys = ene_unphys, ene_phys

    return ene_phys,ene_unphys,min_epsilon,sgn,det_u,det_vt


def main():
    """Output minimum eigenenergy for vortex-free states"""
    global L1, L2, M, Jz
    if len(sys.argv)<4:
        sys.exit("./kitaev.py L1 L2 M [Jz]")
    L1 = int(sys.argv[1])
    L2 = int(sys.argv[2])
    M = int(sys.argv[3])
    Jz = 1.0 if len(sys.argv)<5 else float(sys.argv[4])

    print "# L1={0} L2={1} M={2} Jz={3} Jx=Jy=1.0".format(L1, L2, M, Jz)
    for loop in range(4):
        bond = Bond(loop)
        e = min_energy(bond)[0]
        print "{0}\t{1:.12e}".format(loop, e)


if __name__=="__main__":
    main()
