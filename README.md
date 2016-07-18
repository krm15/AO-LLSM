    # AO-LLSM
Code for supporting AO-LLSM analysis

Compilation: 
A CMakeLists.txt file is provided to compile the project using a C++ compiler. Use CMake and link to ITK libraries. 

Tile stitching using TileReader.cxx.
Name of binary: tileReader
C++ Code stitches across 3D tiles and assembles any 5D region of interest. 
The code is multithreaded and parallelized.

Usage: 
$<path-to-executable>/tileReader iSettingsFile iTileDir oOutputImageDir <OPTIONS> 
Usage: 

iSettings file directory -- Path of CSV file
iTile directories - Path of individual image tiles
oOutputImageDir - Stitched z plane directories 

Options:
-h   --help    Prints this help 
-b   --blend   Off Blends tiles at overlap 
-m   --mip     Off MIP of stitched tiles 
-c   --channel 0   (default) channel value
-t   --time    0   (default) timepoint
-s   --zstart  0   (default) z start plane
-e   --zend    100 (default) z end plane
-n   --threads 1   (default) number of threads
-l   --lsMap       (default) correction filename
-d   --darkLevel 30 (default) correction threshold
-v   --var     2.0 (default) smoothing scale
-x   --exp     _ch (default) string marking channel information
-d   --deconv  ~/  (default) deconvolve tiles based on PSF
-p   --sxy     1   (default) subsampling rate in X/Y
-q   --sz      1   (default) subsampling rate in Z
-o   --offset  ~/  (default) offset filename



Tile registration using TileRegistration.cxx. 
Name of binary: tileRegistration
Usage: 
iSettings file directory -- Path of CSV file
iTile directories - Path of individual image tiles
oOutputImageDir - Stitched z plane directories 

iSettings file directory 
iTile directories 

-h   --help    Prints this help 
-c   --channel 0   (default) channel value
-t   --time    0   (default) timepoint
-n   --threads 1   (default) number of threads
-x   --exp     _ch (default) string marking channel information
-o   --offset  ~/  (default) offset filename
-r   --reg    Off  (default) register z tiles
-s   --zstart  0   (default) z tile


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

