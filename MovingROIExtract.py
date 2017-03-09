#!/usr/bin/python
# Filename : var.py

import sys, os, BatchProcessFunctions

# Default arguments
InputDir = sys.argv[1]
OutputDir = sys.argv[2]
delX = float(sys.argv[3])
delY = float(sys.argv[4])
delZ = float(sys.argv[5])
startX = float(sys.argv[6])
startY = float(sys.argv[7])
startZ = float(sys.argv[8])
sizeX = float(sys.argv[9])
sizeY = float(sys.argv[10])
sizeZ = float(sys.argv[11])
start = int(sys.argv[12])
stop = int(sys.argv[13])
current = os.getcwd()

# Define executable
roiExtract = "/home/ias3/bin/GITROOT/ao-llsm/bin/roiExtract"

for seg in range(start, stop):
    inputImg = InputDir + str(seg) + '_ch0.tif'
    outputImg = OutputDir + str(seg) + '_ch0.tif'

    # Define start
    _startX = startX + delX*seg
    _startY = startY + delY*seg
    _startZ = startZ + delZ*seg
    
    # Define end
    _endX = _startX + sizeX
    _endY = _startY + sizeY
    _endZ = _startZ + sizeZ
    
    # Define command
    Command = "bsub -o log.txt -q megason_2h %s %s %s %s %s %s %s %s %s " % (roiExtract, inputImg, outputImg,
                                               _startX, _startY, _startZ,
                                               _endX, _endY, _endZ)
    print Command
    os.system(Command)
    
    
