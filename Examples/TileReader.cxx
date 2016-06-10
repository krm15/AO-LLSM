/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 667 $  // Revision of last commit
  Date: $Date: 2009-09-16 13:12:21 -0400 (Wed, 16 Sep 2009) $  // Date of last commit
=========================================================================*/

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


int main ( int argc, char* argv[] )
{
  if ( argc < 7 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " iInputSettingsFile iInputImageDir oOutputImageDir ";
    std::cerr << "iChannelNumber iTimePoint iZStart iZEnd" << std::endl;
    return EXIT_FAILURE;
  }

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
  typedef itk::SettingsInfoExtractionFilter< double, ImageType > SettingsFilterType;

  std::string settingsFilename = argv[1];
  std::ifstream infile ( settingsFilename.c_str() );
  if (!infile)
  {
    std::cout << "error in file opening" << std::endl;
    return 0;
  }

  SettingsFilterType::Pointer settingsReader = SettingsFilterType::New();
  settingsReader->SetTileDirectory( argv[2] );
  settingsReader->SetChannelNumber( atoi(argv[4]) );
  settingsReader->SetTimePoint( atoi(argv[5]) );
  settingsReader->Read( infile );
  infile.close();

  StringVectorType m_SettingName;
  DoubleVectorType m_SettingValue;

  settingsReader->GetSettingFieldName( m_SettingName );
  settingsReader->GetSettingFieldValue( m_SettingValue );

  // Setup the dimensions of the largest stitched image
  unsigned int numOfTiles = settingsReader->GetTotalNumberOfTiles();
  unsigned int *tileNumber;
  tileNumber = settingsReader->GetTileNumber();

  double *tileSize;
  tileSize = settingsReader->GetTileSize();

  SizeType tilePixelDimension = settingsReader->GetTileDimension();
  SpacingType spacing = settingsReader->GetTileSpacing();

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


  settingsReader->CreateStitchImage();

  std::cout << std::endl;
  std::cout << "Stitched image origin" << std::endl;
  std::cout << settingsReader->GetStitchOrigin() << std::endl;
  std::cout << "Stitched image dimensions" << std::endl;
  std::cout << settingsReader->GetStitchDimension() << std::endl;
  std::cout << "Stitched image size" << std::endl;
  std::cout << settingsReader->GetStitchSize()[0] << ' ';
  std::cout << settingsReader->GetStitchSize()[1] << ' ';
  std::cout << settingsReader->GetStitchSize()[2] << std::endl;

  // Given zStart and zEnd, assemble an ROI
  unsigned int zStart = atoi( argv[6] );
  unsigned int zEnd = atoi( argv[7] );

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
  ImageType::Pointer m_StitchedImage = settingsReader->GetStitchImage();
  m_StitchedImage->TransformIndexToPhysicalPoint( tempIndex, roiOrigin );

  settingsReader->SetROIOrigin( roiOrigin );
  settingsReader->SetROI( roi );

  std::cout << "Allocate ROI..." << std::endl;

  settingsReader->AllocateROI();

  std::stringstream oFilename;
  oFilename << argv[3] << "ch" << argv[4] << "_" ;
  oFilename << std::setfill( '0' ) << std::setw( 4 ) << argv[5];
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
  series_writer->SetInput( settingsReader->GetROIImage() );
  series_writer->SetFileNames( nameGeneratorOutput->GetFileNames() );
  series_writer->Update();

  return EXIT_SUCCESS;
  }
