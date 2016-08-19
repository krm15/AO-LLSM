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
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include <itkMaximumProjectionImageFilter.h>
#include "itkFillROIImageFilter.h"
#include "itkStitchingSharedData.h"
#include "anyoption.h"
#include "itkTimeProbe.h"

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
  
  
  // Measure time taken
  itk::TimeProbe cputimer;
  cputimer.Start();

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
  opt->addUsage( " -s   --start   0,0,0 (default) start coordinates" );
  opt->addUsage( " -e   --end    100,100,100 (default) end coordinates" );
  opt->addUsage( " -n   --threads 1   (default) number of threads" );
  opt->addUsage( " -l   --lsMap       (default) correction filename" );
  opt->addUsage( " -i   --darkLevel 30 (default) dark level intensity" );
  opt->addUsage( " -v   --var     2.0 (default) smoothing scale" );
  opt->addUsage( " -x   --exp     _ch (default) string marking channel information" );
  opt->addUsage( " -d   --deconv  ~/  (default) deconvolve tiles based on PSF" );
  opt->addUsage( " -p   --sample  1,1 (default) subsampling rate" );
  opt->addUsage( " -o   --offset  ~/  (default) offset filename" );
  opt->addUsage( " -w   --rescale 1.0 (default) rescale factor" );
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
  opt->setOption(  "start",   's' );
  opt->setOption(  "end",     'e' );
  opt->setOption(  "lsMap",   'l' );
  opt->setOption(  "thresh",  'i' );
  opt->setOption(  "var",     'v' );
  opt->setOption(  "threads", 'n' );
  opt->setOption(  "exp",     'x' );
  opt->setOption(  "deconv",  'd' );
  opt->setOption(  "sample",  'p' );
  opt->setOption(  "offset",  'o' );
  opt->setOption(  "rescale", 'w' );

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
  std::vector<unsigned int> start;
  std::vector<unsigned int> end;
  unsigned int zEnd = 10;
  std::string lsMap = "~/";
  std::string PSFImagePath = "~/";
  std::string OffsetFilePath = "~/";
  std::string searchCH = "_ch";
  double thresh = 104.0;
  double var = 100.0;
  unsigned int numOfThreads = 1;
  double rescaleFactor = 1.0;

  std::vector<double> sxy;
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
  if( opt->getValue( 's' ) != NULL  || opt->getValue( "start" ) != NULL  )
  {
    zStart = atoi( opt->getValue( 's' ) );

    std::string input = opt->getValue( 's' );
    std::istringstream ss( input );
    std::string token;

    while( std::getline(ss, token, ',') )
    {
      start.push_back( atoi(token.c_str()) );
    }
  }
  else
  {
    start.push_back( 0 );
    start.push_back( 0 );
    start.push_back( 0 );
  }


  if( opt->getValue( 'e' ) != NULL  || opt->getValue( "end" ) != NULL  )
  {
    zEnd = atoi( opt->getValue( 'e' ) );

    std::string input = opt->getValue( 'e' );
    std::istringstream ss( input );
    std::string token;

    while( std::getline(ss, token, ',') )
    {
      end.push_back( atoi(token.c_str()) );
    }
  }
  else
  {
    end.push_back( 0 );
    end.push_back( 0 );
    end.push_back( 0 );
  }
  if( opt->getValue( 'i' ) != NULL  || opt->getValue( "intensity" ) != NULL  )
  {
    thresh = atof( opt->getValue( 'i' ) );
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
  if( opt->getValue( 'p' ) != NULL  || opt->getValue( "sample" ) != NULL  )
  {
    subsample = true;

    std::string input = opt->getValue( 'p' );
    std::istringstream ss( input );
    std::string token;

    while( std::getline(ss, token, ',') )
    {
      sxy.push_back( atof(token.c_str()) );
    }
  }
  else
  {
    sxy.push_back( 1 );
    sxy.push_back( 1 );
  }

  SharedDataType::Pointer m_SharedData = SharedDataType::New();

  if( opt->getValue( 'w' ) != NULL  || opt->getValue( "rescale" ) != NULL  )
  {
    rescaleFactor = atof( opt->getValue( 'w' ) );
    m_SharedData->m_ScalingFactor = rescaleFactor;
  }
  else
  {
    std::string mipFolder = argv[2];
    std::string searchStringCH = searchCH + opt->getValue( 'c' ) + "_";
    mipFolder += "MIPs/";
    m_SharedData->ComputeScalingFactor( mipFolder, searchStringCH );
  }

  SettingsFilterType::Pointer settingsReader = SettingsFilterType::New();
  settingsReader->SetSettingsDirectory( argv[1] );
  settingsReader->SetTileDirectory( argv[2] );
  settingsReader->SetChannelNumber( ch );
  settingsReader->SetChannelPrefix( searchCH );
  settingsReader->SetTimePoint( tp );
  settingsReader->SetSharedData( m_SharedData );
  settingsReader->Read();

  // Setup the dimensions of the largest stitched image
  unsigned int numOfTiles = settingsReader->GetNumberOfTiles();
  IndexType tileNumber = m_SharedData->m_TileNumber;
  PointType tileSize = m_SharedData->m_TileSize;
  SizeType tilePixelDimension = m_SharedData->m_TileDimension;
  SpacingType spacing = m_SharedData->m_TileSpacing;
  PointType tileOverlap = m_SharedData->m_TileOverlap;

  if( opt->getValue( 'o' ) != NULL  || opt->getValue( "offset" ) != NULL  )
  {
    OffsetFilePath = opt->getValue( 'o' );
    settingsReader->SetOffsetFilePath( OffsetFilePath );
  }
  if( opt->getFlag( "mip" ) || opt->getFlag( 'm' ) )
  {
    mip = true;
  }

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
    std::cout << "Read Correction map with dark level and variance of "  << thresh << ' ' << var << std::endl;
  }

  std::cout << "Sample scan: " << settingsReader->GetSampleScan() << std::endl;
  std::cout << "Scope name : " << settingsReader->GetScopeName() << std::endl;

  std::cout << "Number of tiles " << numOfTiles << std::endl;
  std::cout << "Tile number     " << tileNumber << std::endl;
  std::cout << "Tile size (um)  " << tileSize << std::endl;
  std::cout << "Tile dimension  " << tilePixelDimension << std::endl;
  std::cout << "Tile spacing    " << spacing << std::endl;
  std::cout << "Tile overlap    " << tileOverlap << std::endl;

  ImageType::PointType sOrigin = settingsReader->GetStitchOrigin();

  std::cout << std::endl;
  std::cout << "Stitch origin " << sOrigin << std::endl;
  std::cout << "Stitch dim    " << settingsReader->GetStitchDimension() << std::endl;
  std::cout << "Stitch size   " << settingsReader->GetStitchSize() << std::endl;

  IndexType  roiIndex;
  roiIndex.Fill( 0 );

  IndexType tempIndex;

  SizeType roiSize = settingsReader->GetStitchDimension();
  for( unsigned int i = 0; i < Dimension; i++ )
  {
    if ( end[i] > settingsReader->GetStitchDimension()[i] )
    {
      end[i] = settingsReader->GetStitchDimension()[i] - 1;
    }

    if ( start[i] > end[i] )
    {
      std::cout << "zStart is greater than zEnd" << std::endl;
      return EXIT_SUCCESS;
    }
    roiSize[i] = end[i] - start[i] + 1;

    tempIndex[i] = start[i];
  }

  RegionType roi;
  roi.SetIndex( roiIndex );
  roi.SetSize( roiSize );

  PointType  roiOrigin;
  ImageType::Pointer m_StitchedImage = settingsReader->GetStitchedImage();
  m_StitchedImage->TransformIndexToPhysicalPoint( tempIndex, roiOrigin );

  std::cout << std::endl;
  std::cout << "ROI origin " <<  roiOrigin << std::endl;
  std::cout << "ROI size   " <<  roiSize << std::endl;
  std::cout << "ROI index  " <<  tempIndex << std::endl;

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
  ImageType::SizeType nsize;
  ImageType::PointType nOrigin;
  ImageType::Pointer outputImage = fillROI->GetOutput();
  outputImage->DisconnectPipeline();

  IndexType t_start, t_end;

  for( unsigned int i = 0; i < sxy.size(); i += 2 )
  {
    double sXY = sxy[i];
    double sZ = sxy[i+1];
    nspacing[0] = spacing[0]*sXY;
    nspacing[1] = spacing[1]*sXY;
    nspacing[2] = spacing[2]*sZ;

    for( unsigned int j = 0; j < Dimension; j++ )
    {
      t_start[j] = vcl_ceil( (start[j] * spacing[j])/nspacing[j] );
      t_end[j] = vcl_floor( (end[j] * spacing[j])/nspacing[j] );
      nOrigin[j] = sOrigin[j] + t_start[j] * nspacing[j];
      nsize[j] = t_end[j] - t_start[j] + 1;
    }

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
    resample->SetInput ( outputImage );
    resample->SetDefaultPixelValue ( 0 );
    resample->Update();
    ImageType::Pointer outputImage2 = resample->GetOutput();
    outputImage2->DisconnectPipeline();

    if ( mip )
    {
      MIPFilterType::Pointer mipFilter = MIPFilterType::New();
      mipFilter->SetInput( outputImage2 );
      mipFilter->SetProjectionDimension( 2 );
      mipFilter->SetNumberOfThreads( numOfThreads );
      mipFilter->Update();

      std::stringstream oFilename;
      oFilename << argv[3] << settingsReader->GetChannelName();
      oFilename << "_" << sXY << "_" << sZ << "s";
      oFilename << "_" << ch << "ch" ;
      oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << tp << "t";
      oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << t_start[0];
      oFilename << "-" << std::setfill( '0' ) << std::setw( 4 ) << t_end[0] << "x";
      oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << t_start[1];
      oFilename << "-" << std::setfill( '0' ) << std::setw( 4 ) << t_end[1] << "y";
      oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << t_start[2];
      oFilename << "-" << std::setfill( '0' ) << std::setw( 4 ) << t_end[2];
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
      oFilename << "_" << sXY << "_" << sZ << "s";
      oFilename << "_" << ch << "ch" ;
      oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << tp << "t";
      oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << t_start[0];
      oFilename << "-" << std::setfill( '0' ) << std::setw( 4 ) << t_end[0] << "x";
      oFilename << "_" << std::setfill( '0' ) << std::setw( 4 ) << t_start[1];
      oFilename << "-" << std::setfill( '0' ) << std::setw( 4 ) << t_end[1] << "y";
      oFilename << "_%03dz.tif";

      std::cout << oFilename.str().c_str() << std::endl;

      // Set filename format
      NameGeneratorType::Pointer nameGeneratorOutput = NameGeneratorType::New();
      nameGeneratorOutput->SetSeriesFormat( oFilename.str().c_str() );
      nameGeneratorOutput->SetStartIndex( t_start[2] );
      nameGeneratorOutput->SetEndIndex( t_end[2] );
      nameGeneratorOutput->SetIncrementIndex( 1 );

      // Write out using Series writer
      SeriesWriterType::Pointer series_writer = SeriesWriterType::New();
      series_writer->SetInput( outputImage2 );
      series_writer->SetFileNames( nameGeneratorOutput->GetFileNames() );
      series_writer->SetUseCompression( 0 );
      series_writer->Update();
    }
  }

  cputimer.Stop();
  std::cout << "Stitching took " << cputimer.GetMean() << " seconds" << std::endl;  
  
  //delete opt;

  return EXIT_SUCCESS;
}
