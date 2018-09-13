from __future__ import print_function
import numpy as np
import sys
import time


inFilePref = sys.argv[1]

inFileName = inFilePref + "_EvtNum.txt"
outFileName = inFilePref + "_dupes.txt"


# Load array of (evtNum, file index)
print("Loading", inFileName, "into numpy array...", end='') 
sys.stdout.flush()
stTime = time.clock()
allEvts = np.loadtxt(inFileName, dtype=[('entryIdx', np.uint64), ('fileIdx', np.uint16), \
                                        ('evtNum', np.uint64), ('runNum', np.uint32)])
elTime = time.clock() - stTime
print(elTime, "seconds")


# Get array of duplicates 
print("Processing numpy.unique...", end='')
sys.stdout.flush()
stTime = time.clock()
_, uniqIdx = np.unique(allEvts[['evtNum', 'runNum']], return_index=True)
elTime = time.clock() - stTime
print(elTime, "seconds")
dupEvts = np.delete(allEvts, uniqIdx)


# Print 
print("Found", len(dupEvts), "duplicated events") 
np.savetxt(outFileName, dupEvts, fmt='%i', delimiter='\t')
print("Wrote duplicates to", outFileName)
