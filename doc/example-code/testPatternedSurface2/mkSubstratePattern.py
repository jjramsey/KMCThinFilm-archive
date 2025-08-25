#!/usr/bin/env python
import numpy as np

def makeDomain(E_edge, E_mid, halfNVals):
    
    spacing = 1.0/(halfNVals - 1)

    ramp = []
    for i in range(halfNVals):
        ramp.append(i*spacing)

    nVals = 2*halfNVals;

    domain = np.empty((nVals, nVals))

    preFac = E_mid - E_edge;

    for i in range(nVals):
        for j in range(nVals):
            
            iWrapped = i
            jWrapped = j

            if iWrapped >= halfNVals:
                iWrapped = nVals - 1 - iWrapped

            if jWrapped >= halfNVals:
                jWrapped = nVals - 1 - jWrapped

            rampInd = (iWrapped if (iWrapped < jWrapped) else jWrapped)

            domain[i,j] = preFac*ramp[rampInd] + E_edge

    return domain

# Model A in "COMPUTER SIMULATION OF NUCLEATION ON PATTERNED SURFACES"
# by A. Kuronen, L. Nurminen, and K. Kaski, MRS Fall Meeting 1999
E_edge = 0.65
E_mid = 0.85
halfNVals = 11

domain = makeDomain(E_edge, E_mid, halfNVals)

domainFile = open("singleDomain.dat", 'w')

domainFile.write("%d %d\n" % (domain.shape[0], domain.shape[1]))

for i in range(domain.shape[0]):
    for j in range(domain.shape[1]):
        domainFile.write("%d %d %g\n" % (i, j, domain[i,j]))

domainFile.close()

try:
    from matplotlib.pyplot import imsave
    
    singleSize = 1.0 # inch

    imsave("singleDomain.png", domain, dpi=(domain.shape[1]/singleSize))
except ImportError:
    print("imsave from matplotlib apparently not available. Not saving to PNG.")
