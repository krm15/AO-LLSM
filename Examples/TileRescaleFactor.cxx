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

int main ( int argc, char* argv[] )
{
  const unsigned int Dimension = 3;
  typedef std::vector< std::string > StringVectorType;
  typedef std::vector< double > DoubleVectorType;
  typedef vnl_matrix< double > vnlMatrixType;
  typedef vnl_vector< double > vnlVectorType;

  typedef unsigned short PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;

  typedef ImageType::SpacingType SpacingType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::RegionType RegionType;
  typedef ImageType::PointType PointType;
  typedef SizeType::SizeValueType SizeValueType;

  typedef itk::SettingsInfoExtractionFilter< double, ImageType > SettingsFilterType;
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
  opt->addUsage( " oScaleFactorFile " );
  opt->addUsage( " -h   --help    Prints this help " );
  opt->addUsage( " -c   --channel 0   (default) channel value" );
  opt->addUsage( " -t   --time    0   (default) timepoint" );
  opt->addUsage( " -x   --exp     _ch (default) string marking channel information" );
  opt->addUsage( " -x   --exp     _ch (default) string marking channel information" );
  opt->addUsage( "" );

  /* 4. SET THE OPTION STRINGS/CHARACTERS */

  /* by default all  options will be checked on the command line
    and from option/resource file */

  /* a flag (takes no argument), supporting long and short form */
  opt->setFlag(  "help",  'h' );
  opt->setFlag(  "Physical",  'p' );
  
  /* an option (takes an argument), supporting long and short form */
  opt->setOption(  "channel", 'c' );
  opt->setOption(  "time",    't' );
  opt->setOption(  "exp",     'x' );

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
  std::string searchCH = "_ch";

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
  if( opt->getValue( 'x' ) != NULL  || opt->getValue( "exp" ) != NULL  )
  {
    searchCH = opt->getValue( 'x' );
  }

  SharedDataType::Pointer m_SharedData = SharedDataType::New();

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

  std::cout << "Number of tiles " << numOfTiles << std::endl;
  std::cout << "Tile number     " << tileNumber << std::endl;
  std::cout << "Tile size (um)  " << tileSize << std::endl;
  std::cout << "Tile dimension  " << tilePixelDimension << std::endl;
  std::cout << "Tile spacing    " << spacing << std::endl;

  DirectoryPointer directory = DirectoryType::New();
  directory->Load( argv[2] );

  std::string filename;
  std::stringstream searchString;
  searchString << "_MIP_z.tif";

  unsigned int histogramSize = 500000;
  std::vector< unsigned int > histogram;
  unsigned int m_TrueCountOfTiles = 0;
  for ( unsigned int m = 0; m < directory->GetNumberOfFiles(); m++)
  {
    //std::cout << "m: " << m << std::endl;
    filename = directory->GetFile( m );
    if ( filename.find( searchString.str() ) != std::string::npos )
    {
      std::ifstream infile( filename.c_str() );
      if ( infile )
      {
       infile.close();

       std::cout << filename_mip << std::endl;
       ReaderPointer reader = ReaderType::New();
       reader->SetFileName ( filename.c_str() );

       try
       {
         reader->Update();
       }
       catch( itk::ExceptionObject & err )
       {
         std::cerr << "ExceptionObject caught !" << std::endl;
         std::cerr << err << std::endl;
       }

       ImagePointer img = reader->GetOutput();
       totalPixelCount += img->GetLargestPossibleRegion().GetNumberOfPixels();

       IteratorType It( img, img->GetLargestPossibleRegion() );
       It.GoToBegin();
       while( !It.IsAtEnd() )
       {
         p = static_cast<unsigned int>(It.Get());
         if ( p >= histogramSize )
         {
           p = histogramSize-1;
         }
         histogram[p]++;
         ++It;
      }
      m_TrueCountOfTiles++;
    }
  }
}

  double rescaleFactor = 1.0;
  if ( totalPixelCount > 0 )
  {
    unsigned int topPercentOfPixels = 0.03 * totalPixelCount;
    unsigned int max = histogramSize - 1;
    unsigned int cumsum = 0;
    while( ( max > 0 ) && ( cumsum < topPercentOfPixels ) )
    {
      cumsum += histogram[max];
      max--;
    }

    if ( max > 0 )
    {
     outputFile << double(65535)/( (double)max ) << std::endl;
    }
  }
  std::cout << "Scaling Factor: " << rescaleFactor << std::endl;

  std::ofstream outputFile( argv[2] );
  outputFile << rescaleFactor << std::endl;
  outputFile.close();
  return EXIT_SUCCESS;
}
