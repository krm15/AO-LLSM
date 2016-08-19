#!/usr/bin/python
# Filename : var.py

import sys, os

# Default arguments
TileDir = sys.argv[1]
Prefix = sys.argv[2]
increment = int(sys.argv[3])
start = int(sys.argv[4])
stop = int(sys.argv[5])
channel = sys.argv[6]
current = os.getcwd()

# Define executable
cmd = "/groups/betzig/home/asanos/bin/GITROOT/ao-llsm/release/bin/stitchTiles"

for seg in range(start, stop, increment):
    deconDir =	TileDir + "matlab_decon/"
    stitchDir =	TileDir + "Stitch/" + Prefix
    zStart = seg
    zStop = seg+increment

    # Define command
    Stitch = "qsub -j y -b y -cwd -A hackathon -l sandy=true -V -pe batch 16 %s %s %s %s -t 0 -n 4 -p 1,1 -c %s -s 0,0,%s -e 1000000,1000000,%s" % (cmd, TileDir, deconDir, stitchDir, channel, zStart, zStop)
    
    print Stitch
    os.system(Stitch)
