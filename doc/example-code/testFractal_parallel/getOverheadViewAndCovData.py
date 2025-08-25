#!/usr/bin/env python

from matplotlib.pyplot import imsave
import matplotlib.cm

import glob, os, sys, collections
import numpy as np

myDir = (sys.argv[1] if len(sys.argv) > 1 else '')

if myDir:
    os.chdir(myDir)

inpFileRootPrefix = 'outFile_ProcCoords'
inpFileRootSuffix = '_snapshot'
inpFileSuffix = '.dat'

inpFilePrefix0 = inpFileRootPrefix + '0_0' + inpFileRootSuffix
sliceFor0SnapshotNum = slice(len(inpFilePrefix0), -len(inpFileSuffix))

snapshotNumStrs = []

for file in glob.iglob(inpFilePrefix0 + '*' + inpFileSuffix):
    snapshotNumStrs.append(file[sliceFor0SnapshotNum])

firstSnapShotStr = snapshotNumStrs[0]
sliceForProcCoordStr = slice(len(inpFileRootPrefix), -len(inpFileRootSuffix + firstSnapShotStr + inpFileSuffix))

procCoordStrs = []

istart = []
iendP1 = []
jstart = []
jendP1 = []

for file in glob.iglob(inpFileRootPrefix + '*' + inpFileRootSuffix + firstSnapShotStr + inpFileSuffix):
    procCoordStrs.append(file[sliceForProcCoordStr])

    inpFile = open(file, 'r')
    firstLineFields = inpFile.readline().split()
    inpFile.close()

    istart.append(int(firstLineFields[1]))
    iendP1.append(int(firstLineFields[2]))

    jstart.append(int(firstLineFields[3]))
    jendP1.append(int(firstLineFields[4]))
    
istartGlobal = min(istart)
iendP1Global = max(iendP1)
isize = iendP1Global - istartGlobal

jstartGlobal = min(jstart)
jendP1Global = max(jendP1)
jsize = jendP1Global - jstartGlobal

myCmap = matplotlib.cm.gray

Coverage = collections.namedtuple("Coverage", "cov stdcov snapshot")

myCovData = {}

for snapshotNumStr in snapshotNumStrs:

    img = np.empty((isize, jsize), dtype=np.uint8)

    t = 0.0
    numParticles = 0

    for procCoordStr in procCoordStrs:
        
        fName = inpFileRootPrefix + procCoordStr + inpFileRootSuffix + snapshotNumStr + inpFileSuffix
        
        inpFile = open(fName, 'r')

        for line in inpFile:
            fields = line.split()

            if fields[0] == '#':
                t = float(fields[5].replace('time:', ''))
            else:
                i = int(fields[0])
                j = int(fields[1])

                img[i,j] = np.uint8(fields[2])

                if img[i,j] > 0:
                    numParticles += img[i,j]
        
        inpFile.close()
    
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
