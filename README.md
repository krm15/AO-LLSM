# AO-LLSM
Code for supporting AO-LLSM analysis

Compilation: 
A CMakeLists.txt file is provided to compile the project using a C++ compiler. Use CMake and link to ITK libraries. 

Image stitching using TileReader.cxx. 
C++ Code stitches across 3D tiles and assembles any 5D region of interest. 
The code is multithreaded and parallelized.

Usage: 
$<path-to-executable>/tileReader iInputSettingsFile iInputImageDir oOutputImageDir iChannelNumber iTimePoint iZStart iZEnd 

iInputSettingsFile - CSV File providing information on the acquisition parameters

iInputImageDir - Path to directory containing 3D tiles (any file format) with a specific naming format

oOutputImageDir - Path to directory for outputting individual stitched z-planes

iChannelNumber - Channel number of the tiles to be stitched (0, 1, 2 etc)

iTimePoint - Timpoint at which tiles will be stitched

iZStart - Beginning z-plane for tiling

iZStop - Ending z-plane for tiling



Convert image formats using ConvertFormat.cxx
Usage: 
./convertFormat iInputSettingsFile iInputImageDir iChannelNumber iTimePoint oOutputImageDir


Deconvolve brightfield images using Deconvolution.cxx
Usage: 
./deconvolve iInputImage iKernelImage oOutputImage numOfIterations


Maximum Projection using MaxProjection.cxx
Usage: ./maxProj iInputFile oOutputFile <dim>

Extract region of interest in a timelapse dataset
Usage: 
./roiExtract input output <startX> <startY> <startZ> <endX> <endY> <endZ>

