#!/usr/bin/env python

from matplotlib.pyplot import imsave
import matplotlib.cm

import glob, os, sys, collections
import numpy as np

myDir = (sys.argv[1] if len(sys.argv) > 1 else '')

if myDir:
    os.chdir(myDir)

myCmap = matplotlib.cm.gray

Coverage = collections.namedtuple("Coverage", "cov stdcov snapshot")

myCovData = {}

snapshotStartStr = 'snapshot'
snapshotEndStr = '.dat'
sliceForSnapshotNum = slice(len(snapshotStartStr), -len(snapshotEndStr))

for snapshotFName in glob.iglob(snapshotStartStr + '*' + snapshotEndStr):

    snapshotNumStr = snapshotFName[sliceForSnapshotNum]

    snapshotFile = open(snapshotFName, 'r')

    fields = snapshotFile.readline().split()

    iminGlobal = int(fields[1])
    imaxP1Global = int(fields[2])

    jminGlobal = int(fields[3])
    jmaxP1Global = int(fields[4])

    t = float(fields[5].replace('time:', ''))

    isize = imaxP1Global - iminGlobal
    jsize = jmaxP1Global - jminGlobal

    img = np.empty((isize, jsize), dtype=np.uint8)

    numParticles = 0

    for line in snapshotFile:
        fields = line.split()

        i = int(fields[0])
        j = int(fields[1])

        img[i,j] = np.uint8(fields[2])

        if img[i,j] > 0:
            numParticles += img[i,j]
        
    snapshotFile.close()
    
    imgSz = img.shape[0]*img.shape[1]

    myCovData[t] = Coverage(cov=float(numParticles)/imgSz,
                            stdcov=np.std(img),
                            snapshot=snapshotNumStr)

    imsave('img' + snapshotNumStr.zfill(4) + '.png', img, cmap=myCmap)

covDataFile = open('coverage.dat', 'w')

covDataFile.write("# Snapshot.Num Simulation.Time Coverage Rms.Height\n")

for t in sorted(myCovData.keys()):
    covDataFile.write('%s %g %g %g\n' % (myCovData[t].snapshot, t, myCovData[t].cov, myCovData[t].stdcov))

covDataFile.close()
