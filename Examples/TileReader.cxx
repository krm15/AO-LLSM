/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "itkDirectory.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageSeriesWriter.h"
#include "itkNumericSeriesFileNames.h"

#include "itkSettingsInfoExtractionFilter.h"

#include "itkRichardsonLucyDeconvolutionImageFilter.h"

#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include <itkMaximumProjectionImageFilter.h>
#include "itkFillROIImageFilter.h"
#include "itkStitchingSharedData.h"
#include "anyoption.h"

int main ( int argc, char* argv[] )
{
  const unsigned int Dimension = 3;
  typedef std::vector< std::string > StringVectorType;
  typedef std::vector< double > DoubleVectorType;
  typedef vnl_matrix< double > vnlMatrixType;
  typedef vnl_vector< double > vnlVectorType;
  typedef itk::Directory DirectoryType;
  typedef std::vector< std:: vector< std::vector< std::string > > > StringArray3DType;

  typedef unsigned short PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;

  typedef ImageType::SpacingType SpacingType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::RegionType RegionType;
  typedef ImageType::PointType PointType;
  typedef SizeType::SizeValueType SizeValueType;

  typedef itk::Image< PixelType, 2 > RImageType;
  typedef itk::NumericSeriesFileNames NameGeneratorType;
  typedef itk::ImageSeriesWriter< ImageType, RImageType> SeriesWriterType;
  typedef itk::ImageFileWriter< RImageType> RWriterType;
  typedef itk::SettingsInfoExtractionFilter< double, ImageType > SettingsFilterType;

  typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
  typedef itk::NearestNeighborInterpolateImageFunction< ImageType > InterpolatorType;
  typedef itk::AffineTransform< double, Dimension > TransformType;

  typedef itk::MaximumProjectionImageFilter< ImageType, RImageType > MIPFilterType;

  typedef itk::FillROIImageFilter< ImageType > FillROIFilterType;
  typedef itk::StitchingSharedData< ImageType > SharedDataType;

  /* 1. CREATE AN OBJECT */
  AnyOption *opt = new AnyOption();

  /* 2. SET PREFERENCES  */
  //opt->noPOSIX(); /* do not check for POSIX style character options */
  //opt->setVerbose(); /* print warnings about unknown options */
  //opt->autoUsagePrint(true); /* print usage for bad options */

  /* 3. SET THE USAGE/HELP   */
  opt->addUsage( "" );
  opt->addUsage( "Usage: " );
  opt->addUsage( "" );
  opt->addUsage( " iSettings file directory " );
  opt->addUsage( " iTile directories " );
  opt->addUsage( " oStitched plane directories " );
  opt->addUsage( " -h   --help    Prints this help " );
  opt->addUsage( " -b   --blend   Off Blends tiles at overlap " );
  opt->addUsage( " -m   --mip     Off MIP of stitched tiles " );
  opt->addUsage( " -c   --channel 0   (default) channel value" );
  opt->addUsage( " -t   --time    0   (default) timepoint" );
  opt->addUsage( " -s   --zstart  0   (default) z start plane" );
  opt->addUsage( " -e   --zend    100 (default) z end plane" );
  opt->addUsage( " -n   --threads 1   (default) number of threads" );
  opt->addUsage( " -l   --lsMap       (default) correction filename" );
  opt->addUsage( " -d   --darkLevel 30 (default) correction threshold" );
  opt->addUsage( " -v   --var     2.0 (default) smoothing scale" );
  opt->addUsage( " -x   --exp     _ch (default) string marking channel information" );
  opt->addUsage( " -d   --deconv  ~/  (default) deconvolve tiles based on PSF" );
  opt->addUsage( " -p   --sxy     1   (default) subsampling rate in X/Y" );
  opt->addUsage( " -q   --sz      1   (default) subsampling rate in Z" );
  opt->addUsage( " -o   --offset  ~/  (default) offset filename" );
  opt->addUsage( "" );

  /* 4. SET THE OPTION STRINGS/CHARACTERS */

  /* by default all  options will be checked on the command line
    and from option/resource file */

  /* a flag (takes no argument), supporting long and short form */
  opt->setFlag(  "help",  'h' );
  opt->setFlag(  "blend", 'b' );
  opt->setFlag(  "mip",   'm' );

  /* an option (takes an argument), supporting long and short form */
  opt->setOption(  "channel", 'c' );
  opt->setOption(  "time",    't' );
  opt->setOption(  "zstart",  's' );
  opt->setOption(  "zend",    'e' );
  opt->setOption(  "lsMap",   'l' );
  opt->setOption(  "thresh",  't' );
  opt->setOption(  "var",     'v' );
  opt->setOption(  "threads", 'n' );
  opt->setOption(  "exp",     'x' );
  opt->setOption(  "deconv",  'd' );
  opt->setOption(  "sxy",     'p' );
  opt->setOption(  "sz",      'q' );
  opt->setOption(  "offset",  'o' );

  /* 5. PROCESS THE COMMANDLINE AND RESOURCE FILE */
  /* read options from a  option/resource file with ':'
  separated options or flags, one per line */

  opt->processFile( ".options" );
  opt->processCommandArgs( argc, argv );

  if( ! opt->hasOptions())
  {
    opt->printUsage();
    delete opt;
    return EXIT_FAILURE;
  }

  unsigned int ch = 0;
  unsigned int tp = 0;
  unsigned int zStart = 0;
  unsigned int zEnd = 10;
  std::string lsMap = "~/";
  std::string PSFImagePath = "~/";
  std::string OffsetFilePath = "~/";
  std::string searchCH = "_ch";
  double thresh = 104.0;
  double var = 100.0;
  unsigned int numOfThreads = 1;
  double sxy = 1.0;
  double sz = 1.0;
  bool subsample = false;
  bool mip = false;

  /* 6. GET THE VALUES */
  if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) )
  {
    opt->printUsage();
    delete opt;
    return EXIT_FAILURE;
  }

  if( opt->getArgc() < 3 )
  {
    std::cerr << "Insufficient # of arguments " << opt->getArgc() << std::endl;
    opt->printUsage();
    delete opt;
    return EXIT_FAILURE;
  }

  if( opt->getValue( 'c' ) != NULL  || opt->getValue( "channel" ) != NULL  )
  {
    ch = atoi( opt->getValue( 'c' ) );
  }
  if( opt->getValue( 't' ) != NULL  || opt->getValue( "time" ) != NULL  )
  {
    tp = atoi( opt->getValue( 't' ) );
  }
  if( opt->getValue( 's' ) != NULL  || opt->getValue( "zstart" ) != NULL  )
  {
    zStart = atoi( opt->getValue( 's' ) );
  }
  if( opt->getValue( 'e' ) != NULL  || opt->getValue( "zend" ) != NULL  )
  {
    zEnd = atoi( opt->getValue( 'e' ) );
  }
  if( opt->getValue( 'd' ) != NULL  || opt->getValue( "darkLevel" ) != NULL  )
  {
    thresh = atof( opt->getValue( 'd' ) );
  }
  if( opt->getValue( 'v' ) != NULL  || opt->getValue( "var" ) != NULL  )
  {
    var = atof( opt->getValue( 'v' ) );
  }
  if( opt->getValue( 'n' ) != NULL  || opt->getValue( "threads" ) != NULL  )
  {
    numOfThreads = atof( opt->getValue( 'n' ) );
  }
  if( opt->getValue( 'x' ) != NULL  || opt->getValue( "exp" ) != NULL  )
  {
    searchCH = opt->getValue( 'x' );
  }
  if( opt->getValue( 'p' ) != NULL  || opt->getValue( "sxy" ) != NULL  )
  {
    sxy = atof( opt->getValue( 'p' ) );
    subsample = true;
  }
  if( opt->getValue( 'q' ) != NULL  || opt->getValue( "sz" ) != NULL  )
  {
    sz = atof( opt->getValue( 'q' ) );
    subsample = true;
  }

  SharedDataType::Pointer m_SharedData = SharedDataType::New();

  SettingsFilterType::Pointer settingsReader = SettingsFilterType::New();
  settingsReader->SetSettingsDirectory( argv[1] );
  settingsReader->SetTileDirectory( argv[2] );
  settingsReader->SetChannelNumber( ch );
  settingsReader->SetChannelPrefix( searchCH );
  settingsReader->SetTimePoint( tp );
  settingsReader->SetSharedData( m_SharedData );

  if( opt->getValue( 'o' ) != NULL  || opt->getValue( "offset" ) != NULL  )
  {
    OffsetFilePath = opt->getValue( 'o' );
    settingsReader->SetOffsetFile( OffsetFilePath );
  }
  if( opt->getFlag( "mip" ) || opt->getFlag( 'm' ) )
  {
    mip = true;
  }

  settingsReader->Read();

  StringVectorType m_SettingName = settingsReader->GetSettingFieldName();
  DoubleVectorType m_SettingValue = settingsReader->GetSettingFieldValue();

  // Setup the dimensions of the largest stitched image
  unsigned int numOfTiles = settingsReader->GetNumberOfTiles();
  unsigned int *tileNumber;
  tileNumber = settingsReader->GetTileNumber();

  double *tileSize;
  tileSize = settingsReader->GetTileSize();

  SizeType tilePixelDimension = settingsReader->GetTileDimension();
  SpacingType spacing = settingsReader->GetTileSpacing();

  m_SharedData->SetTileDimension( tilePixelDimension );
  m_SharedData->SetTileSpacing( spacing );

  if( opt->getValue( 'd' ) != NULL  || opt->getValue( "deconv" ) != NULL  )
  {
    PSFImagePath = opt->getValue( 'd' );
    m_SharedData->SetPSFPath( PSFImagePath );
    m_SharedData->SetDeconvolutionIterations( 15 );
    m_SharedData->ReadPSFImage();
    std::cout << "Read PSF image" << std::endl;
  }

  if( opt->getValue( 'l' ) != NULL  || opt->getValue( "lsMap" ) != NULL  )
  {
    lsMap = opt->getValue( 'l' );
    m_SharedData->SetCorrectionInfo( lsMap, var, thresh );
    std::cout << "Read Correction map" << std::endl;
  }

  std::cout << "Number of tiles " << numOfTiles << std::endl;
  std::cout << "Tile number" << std::endl;
  std::cout << tileNumber[0] << ' ' << tileNumber[1] << ' '
                             << tileNumber[2] << std::endl;
  std::cout << " Tile size (um)" << std::endl;
  std::cout << tileSize[0] << ' ' << tileSize[1] << ' '
                           << tileSize[2] << std::endl;
  std::cout << "Tile pixel dimension" << std::endl;
  std::cout << tilePixelDimension << std::endl;
  std::cout << "Tile spacing" << std::endl;
  std::cout << spacing << std::endl;

  //settingsReader->CreateStitchedImage();

  ImageType::PointType sOrigin = settingsReader->GetStitchOrigin();

  std::cout << std::endl;
  std::cout << "Stitched image origin" << std::endl;
  std::cout << sOrigin << std::endl;
  std::cout << "Stitched image dimensions" << std::endl;
  std::cout << settingsReader->GetStitchDimension() << std::endl;
  std::cout << "Stitched image size" << std::endl;
  std::cout << settingsReader->GetStitchSize()[0] << ' ';
  std::cout << settingsReader->GetStitchSize()[1] << ' ';
  std::cout << settingsReader->GetStitchSize()[2] << std::endl;

  if ( zEnd > settingsReader->GetStitchDimension()[2] )
  {
    zEnd = settingsReader->GetStitchDimension()[2] - 1;
  }

  if ( zStart > zEnd )
  {
    std::cout << "zStart is greater than zEnd" << std::endl;
    return EXIT_SUCCESS;
  }

  IndexType  roiIndex;
  roiIndex.Fill( 0 );

  SizeType roiSize = settingsReader->GetStitchDimension();

  roiSize[2] = zEnd - zStart + 1;

  RegionType roi;
  roi.SetIndex( roiIndex );
  roi.SetSize( roiSize );

  IndexType tempIndex;
  tempIndex.Fill( 0 );
  tempIndex[2] = zStart;

  PointType  roiOrigin;
  ImageType::Pointer m_StitchedImage = settingsReader->GetStitchedImage();
  m_StitchedImage->TransformIndexToPhysicalPoint( tempIndex, roiOrigin );

  ImageType::Pointer m_ROIImage = ImageType::New();
  m_ROIImage->SetOrigin( roiOrigin );
  m_ROIImage->SetSpacing( spacing );
  m_ROIImage->SetRegions( roi );
  m_ROIImage->Allocate();
  m_ROIImage->FillBuffer( 0.0 );

  FillROIFilterType::Pointer fillROI = FillROIFilterType::New();
  fillROI->SetInput( m_ROIImage );
  fillROI->SetSharedData( m_SharedData );
  fillROI->InPlaceOn();
  fillROI->SetNumberOfThreads( numOfThreads );

  if( opt->getFlag( "blend" ) || opt->getFlag( 'b' ) )
  {
    //fillROI->SetBlending( 0 );//1
  }
  else
  {
    //fillROI->SetBlending( 0 );
  }

  fillROI->Update();

  ImageType::SpacingType nspacing;
  ImageType::SizeType nsize = roiSize;
  ImageType::PointType nOrigin = roiOrigin;
  ImageType::Pointer outputImage = fillROI->GetOutput();
  if ( subsample )
  {
    nspacing[0] = spacing[0]*sxy;
    nspacing[1] = spacing[1]*sxy;
    nspacing[2] = spacing[2]*sz;

    unsigned int t_zStart = vcl_ceil( (zStart * spacing[2])/nspacing[2] );
    unsigned int t_zEnd = vcl_floor( (zEnd * spacing[2])/nspacing[2] );
    nOrigin[2] = sOrigin[2] + t_zStart * nspacing[2];

    nsize[0] /= sxy;
    nsize[1] /= sxy;
    nsize[2] = t_zEnd - t_zStart + 1;

    // create the resample filter, transform and interpolator
    TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();

    InterpolatorType::Pointer interp = InterpolatorType::New();
    ResampleFilterType::Pointer resample = ResampleFilterType::New();
    resample->SetTransform ( transform );
    resample->SetInterpolator ( interp );
    resample->SetSize ( nsize );
    resample->SetOutputOrigin ( nOrigin );
    resample->SetOutputSpacing ( nspacing );
    resample->SetInput ( fillROI->GetOutput() );
    resample->SetDefaultPixelValue ( 0 );
    resample->Update();
    outputImage = resample->GetOutput();
    outputImage->DisconnectPipeline();

    zStart = t_zStart;
    zEnd = t_zEnd;
  }

  if ( mip )
  {
    MIPFilterType::Pointer mipFilter = MIPFilterType::New();
    mipFilter->SetInput( outputImage );
    mipFilter->SetProjectionDimension( 2 );
    mipFilter->SetNumberOfThreads( numOfThreads );
    mipFilter->Update();

    std::stringstream oFilename;
    oFilename << argv[3] << settingsReader->GetChannelName();
    oFilename << "_" << ch << "ch" ;
    oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << tp << "t";
    oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << zStart;
    oFilename << "-" << std::setfill( '0' ) << std::setw( 4 ) << zEnd;
    oFilename << "z.tif";

    RWriterType::Pointer writer = RWriterType::New();
    writer->SetInput( mipFilter->GetOutput() );
    writer->SetFileName( oFilename.str() );
    writer->Update();
  }
  else
  {
    std::stringstream oFilename;
    oFilename << argv[3] << settingsReader->GetChannelName();
    oFilename << "_" << ch << "ch" ;
    oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << tp << "t";
    oFilename << "_%03dz.tif";

    std::cout << oFilename.str().c_str() << std::endl;

    // Set filename format
    NameGeneratorType::Pointer nameGeneratorOutput = NameGeneratorType::New();
    nameGeneratorOutput->SetSeriesFormat( oFilename.str().c_str() );
    nameGeneratorOutput->SetStartIndex( zStart );
    nameGeneratorOutput->SetEndIndex( zEnd );
    nameGeneratorOutput->SetIncrementIndex( 1 );

    // Write out using Series writer
    SeriesWriterType::Pointer series_writer = SeriesWriterType::New();
    series_writer->SetInput( outputImage );
    series_writer->SetFileNames( nameGeneratorOutput->GetFileNames() );
    series_writer->Update();
  }

  //delete opt;

  return EXIT_SUCCESS;
}
