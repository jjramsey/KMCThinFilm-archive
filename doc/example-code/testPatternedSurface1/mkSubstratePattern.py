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

numDomainsPerEdge = 16

allDomains = np.tile(domain,
                     (numDomainsPerEdge,numDomainsPerEdge))

domainFile = open("tiled%dx%dDomain.dat" % (numDomainsPerEdge,numDomainsPerEdge), 'w')
domainGpFile = open("tiled%dx%dDomain.gpdat" % (numDomainsPerEdge,numDomainsPerEdge), 'w')

domainFile.write("%d %d\n" % (allDomains.shape[0], allDomains.shape[1]))
domainGpFile.write("# i j E_s\n")

for i in range(allDomains.shape[0]):
    for j in range(allDomains.shape[1]):
        domainFile.write("%d %d %g\n" % (i, j, allDomains[i,j]))
        domainGpFile.write("%d %d %g\n" % (i, j, allDomains[i,j]))

    domainGpFile.write("\n")

domainFile.close()

try:
    from matplotlib.pyplot import imsave
    import matplotlib.cm
    
    tiledSize = 3.0 # inch

    imsave("tiled%dx%dDomain.png" % (numDomainsPerEdge,numDomainsPerEdge),
           allDomains, cmap=matplotlib.cm.gray, dpi=(allDomains.shape[1]/tiledSize))
except ImportError:
    print("imsave from matplotlib apparently not available. Not saving to PNG.")
