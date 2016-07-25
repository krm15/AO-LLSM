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
  opt->addUsage( " begin index " );
  opt->addUsage( " end index " );
  opt->addUsage( " -h   --help    Prints this help " );
  opt->addUsage( " -p   --physical Physical or image coordinates " );
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

  StringVectorType m_SettingName = settingsReader->GetSettingFieldName();
  DoubleVectorType m_SettingValue = settingsReader->GetSettingFieldValue();

  // Setup the dimensions of the largest stitched image
  unsigned int numOfTiles = settingsReader->GetNumberOfTiles();
  unsigned int *m_TileNumber;
  m_TileNumber = settingsReader->GetTileNumber();

  double *tileSize;
  tileSize = settingsReader->GetTileSize();

  SizeType tilePixelDimension = settingsReader->GetTileDimension();
  SpacingType spacing = settingsReader->GetTileSpacing();

  m_SharedData->SetTileDimension( tilePixelDimension );
  m_SharedData->SetTileSpacing( spacing );

  std::cout << "Number of tiles " << numOfTiles << std::endl;
  std::cout << "Tile number" << std::endl;
  std::cout << m_TileNumber[0] << ' ' << m_TileNumber[1] << ' '
                             << m_TileNumber[2] << std::endl;
  std::cout << " Tile size (um)" << std::endl;
  std::cout << tileSize[0] << ' ' << tileSize[1] << ' '
                           << tileSize[2] << std::endl;
  std::cout << "Tile pixel dimension" << std::endl;
  std::cout << tilePixelDimension << std::endl;
  std::cout << "Tile spacing" << std::endl;
  std::cout << spacing << std::endl;

  // beginCorner or endCorner
  IndexType beginIndex, endIndex;
  PointType beginCorner, endCorner;
 
  if( opt->getFlag( "physical" ) || opt->getFlag( 'p' ) )
  {
    beginCorner[0] = atoi( argv[3] );
    beginCorner[1] = atoi( argv[4] );
    beginCorner[2] = atoi( argv[5] );
    endCorner[0] = atoi( argv[6] );
    endCorner[1] = atoi( argv[7] );
    endCorner[2] = atoi( argv[8] );
    
    settingsReader->GetStitchedImage()->TransformPhysicalPointToIndex( beginCorner, beginIndex );
    settingsReader->GetStitchedImage()->TransformPhysicalPointToIndex( endCorner, endIndex );
  }
  else
  {
    beginIndex[0] = atoi( argv[3] );
    beginIndex[1] = atoi( argv[4] );
    beginIndex[2] = atoi( argv[5] );
    endIndex[0] = atoi( argv[6] );
    endIndex[1] = atoi( argv[7] );
    endIndex[2] = atoi( argv[8] );
  
    settingsReader->GetStitchedImage()->TransformIndexToPhysicalPoint( beginIndex, beginCorner );
    settingsReader->GetStitchedImage()->TransformIndexToPhysicalPoint( endIndex, endCorner );
  }

  std::cout << "Begin corner: " << beginCorner << std::endl;
  std::cout << "End corner: " << beginCorner << std::endl;
  
  IndexType m_ScanStart, m_ScanEnd;
  
  std::cout << "Setting scan start and end values for ROI" << std::endl;
  for( unsigned int k = 0; k < Dimension; k++ )
  {
    //std::cout <<  beginCorner << ' ' << endCorner << std::endl;
    double scanStartVal = 100000, scanEndVal = -100000;
    for( unsigned int i = 0; i < m_TileNumber[k]; i++ )
    {
      //std::cout << k << ' ' << i << ' ' << m_SharedData->m_TileCoverStart[k][i] << ' '
      //          << m_SharedData->m_TileCoverEnd[k][i] << std::endl;
      if ( ( beginCorner[k] >= m_SharedData->m_TileCover[k][0][0][i] - 0.0001 ) &&
           ( beginCorner[k] <= m_SharedData->m_TileCover[k][1][0][i] + 0.0001 ) &&
           ( scanStartVal >= m_SharedData->m_TileCover[k][0][0][i] ) )
      {
        m_ScanStart[k] = i;
        scanStartVal =  m_SharedData->m_TileCover[k][0][0][i];
      }

      if ( ( endCorner[k] >= m_SharedData->m_TileCover[k][0][0][i] - 0.001 ) &&
           ( endCorner[k] <= m_SharedData->m_TileCover[k][1][0][i] + 0.001 ) &&
           ( scanEndVal <= m_SharedData->m_TileCover[k][1][0][i] ) )
      {
        m_ScanEnd[k] = i;
        scanEndVal =  m_SharedData->m_TileCover[k][1][0][i];
      }
    }

    unsigned int temp;
    if ( m_ScanStart[k] > m_ScanEnd[k] )
    {
      temp = m_ScanEnd[k];
      m_ScanEnd[k] = m_ScanStart[k];
      m_ScanStart[k] = temp;
    }
  }
  
  std::cout << "X tiles: " << m_ScanStart[0] << ' ' << m_ScanEnd[0] << std::endl;
  std::cout << "Y tiles: " << m_ScanStart[1] << ' ' << m_ScanEnd[1] << std::endl;
  std::cout << "Z tiles: " << m_ScanStart[2] << ' ' << m_ScanEnd[2] << std::endl;
      
  
  for( unsigned int i = m_ScanStart[0]; i <= m_ScanEnd[0]; i++ )
  {
    for( unsigned int j = m_ScanStart[1]; j <= m_ScanEnd[1]; j++ )
    {
      for( unsigned int k = m_ScanStart[2]; k <= m_ScanEnd[2]; k++ )
      {
        //std::cout << counter << ' ' << counter%(m_NumOfValidThreads) << std::endl;
        std::cout << m_SharedData->m_TileFileNameArray[i][j][k] << std::endl;  
      }
    }
  }
  //delete opt;

  return EXIT_SUCCESS;
}
