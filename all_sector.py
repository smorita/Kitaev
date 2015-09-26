#!/usr/bin/env python

"""Search the ground state of the Kitaev model."""

__author__ = "Satoshi Morita"
__version__ = "2.0.0"

import sys
import logging
from itertools import product,combinations
import numpy as np
import kitaev

def shift(vortex,dr):
    from kitaev import position,index
    result = [ 0 for i in range(len(vortex)) ]
    for i,v in enumerate(vortex):
        x,y = position(i)
        result[index(x,y,dr)] = v
    return tuple(result)


def is_valid(vortex):
    total = sum(vortex)
    if total%2 != 0: return False
    elif total == 0: return True
    elif total == len(vortex): return True
    elif vortex[0] > 0: return False
    for i in [ i for i,v in enumerate(vortex) if v==0 ]:
        if i==0: continue
        dr = [ -cood for cood in kitaev.position(i)]
        if vortex>shift(vortex,dr): return False
    return True


def flip(bond,vortex):
    r = [ i for i,v in enumerate(vortex) if v==1 ]
    for i in range(0,len(r),2):
        bond.create_vortex(r[i],r[i+1])


def main():
    """Output the ground state"""
    if len(sys.argv)<4:
        sys.exit("Usage: ./script.py L1 L2 M [Jz]")
    kitaev.L1 = int(sys.argv[1])
    kitaev.L2 = int(sys.argv[2])
    kitaev.M = int(sys.argv[3])
    kitaev.Jz = 1.0 if len(sys.argv)<5 else float(sys.argv[4])

    N_unit = kitaev.L1*kitaev.L2
    logging.basicConfig(format="%(levelname)s:\t%(message)s", level=logging.INFO)

    min_e = (0,0,0)
    for vortex in product((0,1),repeat=N_unit):
        if not is_valid(vortex): continue

        for loop in range(4):
            bond = kitaev.Bond(loop)
            flip(bond,vortex)
            e = kitaev.min_energy(bond)[0]
            logging.info("{0:.12e} {1} {2}".format(e, loop, vortex))
            if min_e[0]>e:
                min_e = (e, loop, vortex)

    output = []
    output.append("########################################")
    output.append("# L1\t=\t{0}".format(kitaev.L1))
    output.append("# L2\t=\t{0}".format(kitaev.L2))
    output.append("# M\t=\t{0}".format(kitaev.M))
    output.append("# Jx\t=\t{0}".format(kitaev.Jx))
    output.append("# Jy\t=\t{0}".format(kitaev.Jy))
    output.append("# Jz\t=\t{0}".format(kitaev.Jz))
    output.append("########################################")
    output.append("# loop\t=\t{0}".format(min_e[1]))
    output.append("# vortex\t=\t{0}".format(min_e[2]))
    output.append("########################################")
    output.append("# Ground-state energy")
    output.append("{0:.12e}".format(min_e[0]))
    print "\n".join(output)


if __name__=="__main__":
    main()
