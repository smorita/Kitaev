#!/usr/bin/env python
import sys
import numpy as np

def index(r, dr=(0,0)):
    x0 = r[0]+dr[0]
    y0 = r[1]+dr[1]
    x = (x0-M*(y0/L2))%L1
    y = y0%L2
    return x+L1*y

def position(i, dr=(0,0)):
    if i<0 or i>=Nunit:
        sys.stderr.write("position(): bad index {0}\n".format(i))
        return
    x0 = (i%L1)+dr[0]
    y0 = (i/L1)+dr[1]
    x = (x0-M*(y0/L2))%L1
    y = y0%L2
    return (x, y)

def keyBondZ(r):
    ia = ib = index(r)
    return (ia,ib,'z')

def keyBondX(r):
    ia = index(r)
    ib = index(r, (1,0))
    return (ia,ib,'x')

def keyBondY(r):
    ia = index(r)
    ib = index(r, (0,1))
    return (ia,ib,'y')

def initBond(loop=0):
    """Create bond"""
    bond = {}
    ## z bond
    for r in np.ndindex(L1,L2):
        bond[ keyBondZ(r) ] = Jz
    ## x bond
    for r in np.ndindex(L1,L2):
        bond[ keyBondX(r) ] = Jx
    ## y bond
    for r in np.ndindex(L1,L2):
        bond[ keyBondY(r) ] = Jy

    if loop==0: pass
    elif loop==1: flipLoop1(bond)
    elif loop==2: flipLoop2(bond)
    elif loop==3: 
        flipLoop1(bond)
        flipLoop2(bond)
    else:
        sys.stderr.write("initBond(): bad loop idx {0}\n".format(loop))
    return bond

def flipLoop1(bond):
    """Insert a Wilson loop along L1*e1"""
    # flip y-bond
    for x in range(L1):
        bond[ keyBondY((x,0)) ] *= -1

def flipLoop2(bond):
    """Insert a Wilson loop along M*e1+L2*e2"""
    # flip x-bond
    for y in range(L2):
        bond[ keyBondX((0,y)) ] *= -1
    # flip y-bond
    for x in range(M):
        bond[ keyBondY((x+1,L2-1)) ] *= -1

def createVortex1(bond,x1,x2,y):
    """create vortices at (x1,y) and (x2,y)"""
    if x1>x2: x1,x2 = x2,x1
    for x in range(x1,x2):
        bond[ keyBondY((x+1,y)) ] *= -1

def createVortex2(bond,x,y1,y2):
    """create vortices at (x,y1) and (x,y2)"""
    if y1>y2: y1,y2 = y2,y1
    for y in range(y1,y2):
        bond[ keyBondX((x,y+1)) ] *= -1

def createVortex(bond,i,j):
    """create vortices at hexagones i and j."""
    xi,yi = position(i)
    xj,yj = position(j)
    createVortex1(bond,xi,xj,yi) # from (xi,yi) to (xj,yi)
    createVortex2(bond,xj,yi,yj) # from (xj,yi) to (xj,yj)

def createVortexConfig(bond,vortexConfig):
    """flip bond such that it has the given vortex configuration"""
    if sum(vortexConfig)%2 != 0:
        sys.stderr.write("createVortexConfig(): bad vortex config "+stringVortex(vortexConfig)+"\n")
    hexagone = []
    for i,v in enumerate(vortexConfig):
        if v%2==1: hexagone.append(i)
    for i in range(0,len(hexagone),2):
        createVortex(bond,hexagone[i],hexagone[i+1])

def stringVortex(vortexConfig):
    string=""
    for v in vortexConfig:
        string += "{0} ".format(v)
    return string

def minEnergy(bond):
    """Calculate minimum energy."""
    EPS = 1.0e-6 # threashold for zero eigenvalue

    # create matrix A
    matA = np.zeros((Nunit,Nunit),dtype=float)
    for (r,v) in bond.iteritems():
        ia,ib,direction = r
        matA[ia][ib] += v
    # B = A^T*A
    matB = np.dot(matA.transpose(), matA)
    # B2 = A*A^T
    matB2 = np.dot(matA, matA.transpose())
    # calculate eigenvalues and eigenvectors
    val,vec = np.linalg.eigh(matB)
    val2,vec2 = np.linalg.eigh(matB2)
    
    # one-particle energy
    epsilon = np.sqrt(np.abs(val))
    
    # energy of the vacuum state
    sumE = 0.0;
    for i,e in enumerate(epsilon):
        if e>EPS: sumE += e
        else: epsilon[i] = 0.0
    sumE *= -0.5

    # create orthogonal transition matricies
    matU = vec
    matV = np.zeros_like(matU)
    for i,e in enumerate(epsilon):
        if e>EPS:
            matV[:,i]=np.dot(matA,vec[:,i])/e
        else:
            matV[:,i]=vec2[:,i]
    # their determinant
    detU = np.linalg.det(matU)
    detV = np.linalg.det(matV)
    
    # calculate parity of the projection operator
    ## product of u_{ij}
    sgn = +1
    for b in bond.itervalues():
        if b<0: sgn *= -1
    ## from boundary condition
    if (L1+L2+M*(L1-M))%2 != 0: sgn *= -1 # (-1)^theta
    ## det(Q) = det(VU)
    if detV*detU < 0: sgn *= -1
    
    # judge whether the vacuume state is physical or not
    minEpsilon = min(epsilon)
    if sgn > 0: # The vacuum state is physical.
        physE = sumE
        unphysE = sumE + minEpsilon
    else: # The vacuum state is unphysical.
        physE = sumE + minEpsilon
        unphysE = sumE

    return physE,unphysE,sgn,detV,detU,minEpsilon

def energy(l1,l2,m,jz=1.0,loop=0,vortex=[]):
    """Calculate lowest energy with given boundary condition and vortex configuration"""
    global L1,L2,M,Jz,Nunit
    L1,L2,M,Jz = l1,l2,m%L1,jz
    Nunit = L1*L2
    bond = initBond(loop)
    if len(vortex)==Nunit:
        createVortexConfig(bond,vortex)
    return minEnergy(bond)

L1 = 6
L2 = 2
M  = 2
Jx = Jy = Jz = 1.0
Nunit = L1*L2

if __name__ == "__main__":
    if len(sys.argv)<4:
        sys.stderr.write("./script.py L1 L2 M [Jz=1.0] \n")
        sys.exit(1)
    
    # L1, L2, M = 6, 2, 2 # hexagon
    # L1, L2, M = 12, 1, 7 # 4x3 rectangle
    L1 = int(sys.argv[1])
    L2 = int(sys.argv[2])
    M  = int(sys.argv[3])
    
    if M >= L1: 
        sys.stderr.write("error: M should be smaller than L1. (L1={0} L2={1} M={2})\n".format(L1,L2,M))
        sys.exit(1)
    
    Jx = Jy = 1.0
    Jz = float(sys.argv[4]) if len(sys.argv)>4 else 1.0

    print "# L1={0} L2={1} M={2} Jz={3} Jx=Jy=1.0".format(L1,L2,M,Jz)
    for l in range(4):
        result = energy(L1,L2,M,Jz,loop=l)
        print  "{0}\t{1:+.12e}".format(l,result[0])
