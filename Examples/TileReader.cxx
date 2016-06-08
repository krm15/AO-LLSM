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


int main ( int argc, char* argv[] )
{
  if ( argc < 6 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " iInputSettingsFile iInputImageDir iChannelNumber ";
    std::cerr << "iTimePoint iZStart iZEnd oOutputImageDir" << std::endl;
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

  typedef ImageType::SizeType SizeType;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::RegionType RegionType;
  typedef ImageType::PointType PointType;
  typedef SizeType::SizeValueType SizeValueType;

  typedef itk::Image< PixelType, 2 > RImageType;
  typedef itk::NumericSeriesFileNames NameGeneratorType;
  typedef itk::ImageSeriesWriter< ImageType, RImageType> SeriesWriterType;

  std::string settingsFilename = argv[1];
  std::ifstream infile ( settingsFilename.c_str() );
  if (!infile)
  {
    std::cout << "error in file opening" << std::endl;
    return 0;
  }

  std::string value, line;
  StringVectorType m_SettingName;
  m_SettingName.resize( 100 );

  DoubleVectorType m_SettingValue;
  m_SettingValue.resize( 100 );

  // Read first two lines
  std::getline ( infile, line);
  std::stringstream nameStream(line);

  std::getline ( infile, line);
  std::stringstream valueStream( line );

  // First two lines is 29 fields of data
  for( unsigned int i = 0; i < 29; i++ )
  {
    std::getline ( nameStream, value, ',' );
    m_SettingName[i] = value;
    //std::cout << value << ' ';
    std::getline ( valueStream, value, ',' );
    m_SettingValue[i] = atof( value.c_str() );
    //std::cout << value << std::endl;
  }

  // Setup the dimensions of the largest stitched image
  unsigned int numOfTiles = 1;
  double tileNumber[3];
  double tileSize[3];

  for( unsigned int i = 0; i < Dimension; i++ )
  {
    numOfTiles *= m_SettingValue[i];
    tileNumber[i] = m_SettingValue[i];
    tileSize[i] = m_SettingValue[9+i];
  }

  // Read next two lines
  StringVectorType m_TileInfoName;
  m_TileInfoName.resize( 100 );

  std::getline ( infile, line);
  std::stringstream tileInfoNameStream( line );

  for( unsigned int i = 0; i < 9; i++ )
  {
    std::getline ( tileInfoNameStream, value, ',' );
    m_TileInfoName[i] = value;
    //std::cout << value << std::endl;
  }

  vnlMatrixType m_TileInfoValue (numOfTiles, 9);
  for( unsigned int i = 0; i < numOfTiles; )
  {
    //    std::cout << i << std::endl;
    std::getline ( infile, line);

    char dummy = line.c_str()[0];
    if( ( dummy != '-' ) )
    {
      std::stringstream tileInfoValueStream( line );

      for( unsigned int j = 0; j < 9; j++ )
      {
        std::getline ( tileInfoValueStream, value, ',' );
        m_TileInfoValue[i][j] = atof( value.c_str() );
        //std::cout << ' ' << value;
      }

      //std::cout << std::endl;
      i++;
    }
    else
    {
      //std::cout << std::endl;
    }
  }

  vnlVectorType tileCenter, stitchCenter, newCenter;
  tileCenter.set_size( Dimension );
  stitchCenter.set_size( Dimension );
  for( unsigned int i = 0; i < Dimension; i++ )
  {
    stitchCenter[i] = m_TileInfoValue.get_column(i+6).mean();
    //std::cout << stitchCenter[i] << ' ';
  }
  //std::cout << std::endl;

  double theta = 31.8 * vnl_math::pi_over_180;
  vnlMatrixType rotMatrix( Dimension, Dimension );
  rotMatrix[0][0] = rotMatrix[0][2] = rotMatrix[1][1] = rotMatrix[2][1] = 0.0;
  rotMatrix[0][1] = -1;
  rotMatrix[1][0] = -vcl_cos( theta );
  rotMatrix[2][2] = vcl_cos( theta );
  rotMatrix[1][2] = rotMatrix[2][0] = vcl_sin( theta );

  //std::cout << rotMatrix << std::endl;

  vnlMatrixType m_TransformedTileInfoValue (numOfTiles, 3*Dimension);
  for( unsigned int i = 0; i < numOfTiles; i++ )
  {
    for( unsigned int j = 0; j < Dimension; j++ )
    {
      tileCenter[j] = m_TileInfoValue[i][j+6];
      //std::cout << m_TileInfoValue[i][j] << ' ';
    }

    newCenter = rotMatrix * ( tileCenter - stitchCenter );

    //std::cout << newCenter << std::endl;

    for( unsigned int j = 0; j < Dimension; j++ )
    {
      m_TransformedTileInfoValue[i][j] = newCenter[j];
      m_TransformedTileInfoValue[i][j+Dimension] = newCenter[j] - 0.5*tileSize[j];
      m_TransformedTileInfoValue[i][j+2*Dimension] = newCenter[j] + 0.5*tileSize[j];
    }
  }

  // Create a vector of tile origins along each axis
  DoubleVectorType tileCoverStart[3];
  DoubleVectorType tileCoverEnd[3];
  for( unsigned int i = 0; i < Dimension; i++ )
  {
    tileCoverStart[i].resize( tileNumber[i] );
    tileCoverEnd[i].resize( tileNumber[i] );
  }

  for( unsigned int i = 0; i < numOfTiles; i++ )
  {
    for( unsigned int j = 0; j < Dimension; j++ )
    {
      unsigned int temp =  m_TileInfoValue[i][j];
      tileCoverStart[j][temp] = m_TransformedTileInfoValue[i][j+3];
      tileCoverEnd[j][temp] = m_TransformedTileInfoValue[i][j+6];
    }
  }

  vnlVectorType minStart( 3 );
  vnlVectorType maxEnd( 3 );
  for( unsigned int i = 0; i < Dimension; i++ )
  {
    minStart[i] = m_TransformedTileInfoValue.get_column(i+3).min_value();
    maxEnd[i] = m_TransformedTileInfoValue.get_column(i+6).max_value();
  }

  infile.close();

  // Read all the files in the input directory of type ch and at timePoint
  std::string filename;
  std::stringstream searchStringCH, searchStringXYZT, filename2;

  DirectoryType::Pointer directory = DirectoryType::New();
  directory->Load( argv[2] );

  StringArray3DType tileFileNameArray;
  tileFileNameArray.resize( m_SettingValue[0] );

  for( unsigned int i = 0; i < tileNumber[0]; i++ )
  {
    tileFileNameArray[i].resize( tileNumber[1] );
    for( unsigned int j = 0; j < tileNumber[1]; j++ )
    {
      tileFileNameArray[i][j].resize( tileNumber[2] );
      for( unsigned int k = 0; k < tileNumber[2]; k++ )
      {
        tileFileNameArray[i][j][k] = std::string();
        searchStringCH << "_ch" << argv[3];
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << i << "x_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << j << "y_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << k << "z_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 4 ) << argv[4] << "t.mha";

        //std::cout << i << j << k << std::endl;
        for ( unsigned m = 0; m < directory->GetNumberOfFiles(); m++)
        {
          filename = directory->GetFile( m );

          if ( ( filename.find( searchStringCH.str() ) != std::string::npos ) &&
               ( filename.find( searchStringXYZT.str() ) != std::string::npos ) )
          {
            //std::cout << filename << std::endl;
            filename2 << argv[2] << filename;
            tileFileNameArray[i][j][k] = filename2.str();
          }
          filename2.str( std::string() );
        }
        searchStringCH.str( std::string() );
        searchStringXYZT.str( std::string() );
      }
    }
  }

  // Read one image to get tilePixelDimensions and spacing
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( tileFileNameArray[0][0][0] );
  reader->SetGlobalWarningDisplay( 0 );
  reader->Update();
  ImageType::Pointer currentImage = reader->GetOutput();
  ImageType::SizeType tilePixelDimension = currentImage->GetLargestPossibleRegion().GetSize();

  ImageType::SpacingType spacing;
  spacing[0] = tileSize[0]/tilePixelDimension[0];
  spacing[1] = tileSize[1]/tilePixelDimension[1];
  spacing[2] = tileSize[2]/tilePixelDimension[2];

  std::cout << "Tile number" << std::endl;
  std::cout << tileNumber[0] << ' ' << tileNumber[1] << ' ' << tileNumber[2] << std::endl;
  std::cout << " Tile size (um)" << std::endl;
  std::cout << tileSize[0] << ' ' << tileSize[1] << ' ' << tileSize[2] << std::endl;
  std::cout << "Tile pixel dimensions" << std::endl;
  std::cout << tilePixelDimension << std::endl;
  std::cout << "Tile spacing" << std::endl;
  std::cout << spacing << std::endl;

  // Create the dimensions of the large image
  double                stitchSize[3];
  ImageType::PointType  stitchOrigin;
  ImageType::SizeType   stitchDimension;
  ImageType::IndexType  stitchIndex;
  ImageType::RegionType stitchRegion;

  for( unsigned int i = 0; i < Dimension; i++ )
  {
    stitchIndex[i]     = 0;
    stitchSize[i]      = maxEnd[i] - minStart[i];
    stitchOrigin[i]    = minStart[i];
    stitchDimension[i] = stitchSize[i]/spacing[i];
  }

  ImageType::Pointer stitchedImage = ImageType::New();
  stitchedImage->SetOrigin( stitchOrigin );
  stitchedImage->SetSpacing( spacing );
  stitchedImage->SetRegions( stitchRegion );

  std::cout << std::endl;
  std::cout << "Stitched image origin" << std::endl;
  std::cout << stitchOrigin << std::endl;
  std::cout << "Stitched image dimensions" << std::endl;
  std::cout << stitchDimension << std::endl;
  std::cout << "Stitched image size" << std::endl;
  std::cout << stitchSize[0] << ' ' << stitchSize[1] << ' ' << stitchSize[2] << std::endl;

  // Given zStart and zEnd, assemble an ROI
  unsigned int zStart = atoi( argv[5] );
  unsigned int zEnd = atoi( argv[6] );

  ImageType::RegionType roi;

  ImageType::IndexType  roiIndex, tempIndex;
  roiIndex = stitchIndex;

  tempIndex = stitchIndex;
  tempIndex[2] = zStart;

  ImageType::SizeType   roiSize;
  roiSize = stitchDimension;
  roiSize[2] = zEnd - zStart + 1;

  roi.SetIndex( roiIndex );
  roi.SetSize( roiSize );

  ImageType::PointType  roiOrigin;
  stitchedImage->TransformIndexToPhysicalPoint( tempIndex, roiOrigin );

  ImageType::Pointer roiImage = ImageType::New();
  roiImage->SetOrigin( roiOrigin );
  roiImage->SetSpacing( spacing );
  roiImage->SetRegions( roi );
  roiImage->Allocate();
  roiImage->FillBuffer( 0.0 );

  // Identify all the tiles that belong to this roi
  double zBeginOrigin = roiOrigin[2];
  double zEndOrigin = roiOrigin[2] + roiSize[2]*spacing[2];

  unsigned int zScanStart = 1000000, zScanEnd = 0;
  for( unsigned int i = 0; i < tileNumber[2]; i++ )
  {
    //std::cout << tileCoverStart[2][i] << ' ' << tileCoverEnd[2][i] << std::endl;

    if ( ( zBeginOrigin >= tileCoverStart[2][i] ) &&
         ( zBeginOrigin <= tileCoverEnd[2][i] ) &&
         ( zScanStart > i ) )
    {
      zScanStart = i;
    }

    if ( ( zEndOrigin >= tileCoverStart[2][i] ) &&
         ( zEndOrigin <= tileCoverEnd[2][i] ) &&
         ( zScanEnd < i ) )
    {
      zScanEnd = i;
    }
  }

  std::cout << std::endl;
  std::cout << "Tile z steps under consideration" << std::endl;
  std::cout << zScanStart << ' ' << zScanEnd << std::endl;

  // Start a loop that will read all the tiles from zScanStart to zScanEnd
  ImageType::PointType currentTileOrigin;
  ImageType::RegionType currentTileRegion, roiSubRegion;
  for( unsigned int i = 0; i < tileNumber[0]; i++ )
  {
    currentTileOrigin[0] = tileCoverStart[0][i];
    for( unsigned int j = 0; j < tileNumber[1]; j++ )
    {
      currentTileOrigin[1] = tileCoverStart[1][j];
      for( unsigned int k = zScanStart; k < zScanEnd; k++ )
      {
        currentTileOrigin[2] = tileCoverStart[2][k];

        filename = tileFileNameArray[i][j][k];
//        std::cout << filename << std::endl;

        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName( filename.c_str() );
        reader->SetGlobalWarningDisplay( 0 );
        reader->Update();

        ImageType::Pointer currentImage = reader->GetOutput();

 //       std::cout << "Current Tile Origin" << std::endl;
 //       std::cout << currentTileOrigin << std::endl;
 //       std::cout << currentImage->GetLargestPossibleRegion().GetIndex() << std::endl;

 //       std::cout << "ROI Origin" << std::endl;
 //       std::cout << roiOrigin << std::endl;
 //       std::cout << roiIndex << std::endl;

        currentImage->SetOrigin( currentTileOrigin );
        currentTileRegion = currentImage->GetLargestPossibleRegion();

        // Map the two regions roiSubRegion, roiOrigin
        SizeType sizeA, sizeB, s;
        sizeA = currentTileRegion.GetSize();
        sizeB = roiSize;

        IndexType tIndexA, tIndexB;
        currentImage->TransformPhysicalPointToIndex( roiOrigin, tIndexA );
        roiImage->TransformPhysicalPointToIndex( currentTileOrigin, tIndexB );

//        std::cout << "Current Tile Index" << std::endl;
//        std::cout << tIndexA << std::endl;

//        std::cout << "ROI Index" << std::endl;
//        std::cout << tIndexB << std::endl;


        IndexType sIndexA, sIndexB;
        for( unsigned int p = 0; p < Dimension; p++ )
        {
          if ( currentTileOrigin[p] > roiOrigin[p] )
          {
            sIndexA[p] = 0;
            sIndexB[p] = tIndexB[p];
            s[p] = sizeA[p];
            if ( s[p] > static_cast< SizeValueType >( sizeB[p] - sIndexB[p] - 1 ) )
            {
              s[p] = sizeB[p] - sIndexB[p];
            }
          }
          else
          {
            sIndexB[p] = 0;
            sIndexA[p] = tIndexA[p];
            s[p] = sizeB[p];
            if ( s[p] > static_cast< SizeValueType >( sizeA[p] - sIndexA[p] - 1 ) )
            {
              s[p] = sizeA[p] - sIndexA[p];
            }
          }
//          std::cout << "Current Tile Index" << std::endl;
//          std::cout << sIndexA << std::endl;

//          std::cout << "ROI Index" << std::endl;
//          std::cout << sIndexB << std::endl;
        }

        currentTileRegion.SetIndex( sIndexA );
        currentTileRegion.SetSize( s );
        roiSubRegion.SetIndex( sIndexB );
        roiSubRegion.SetSize( s );

//        std::cout << roiSubRegion << std::endl;
//        std::cout << currentTileRegion << std::endl;

        // Using these images, fill up roiImage
        IteratorType rIt( roiImage, roiSubRegion );
        rIt.GoToBegin();

        IteratorType cIt( currentImage, currentTileRegion );
        cIt.GoToBegin();

        while( !cIt.IsAtEnd() )
        {
          rIt.Set( cIt.Get() );
          ++cIt;
          ++rIt;
        }
      }
    }
  }

  std::stringstream oFilename;
  oFilename << argv[7] << "ch" << argv[3] << "_" ;
  oFilename << std::setfill( '0' ) << std::setw( 4 ) << argv[4];
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
  series_writer->SetInput( roiImage );
  series_writer->SetFileNames( nameGeneratorOutput->GetFileNames() );
  series_writer->Update();



  return EXIT_SUCCESS;
  }
