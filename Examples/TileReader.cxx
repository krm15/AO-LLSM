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

int main ( int argc, char* argv[] )
  {
  if ( argc < 6 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputSettingsFile inputImageDir channelNumber timePoint output" << std::endl;
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

  std::string settingsFilename = argv[1];
  std::ifstream infile ( settingsFilename );
  if (!infile)
  {
    std::cout << "error in file opening" << std::endl;
    return 0;
  }

  std::string value, line;
  StringVectorType m_SettingName;
  m_SettingName.reserve( 100 );

  DoubleVectorType m_SettingValue;
  m_SettingValue.reserve( 100 );

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
       //std::cout << value << std::endl;
       std::getline ( valueStream, value, ',' );
       m_SettingValue[i] = atof( value.c_str() );
       //std::cout << value << std::endl;
  }

  // Setup the dimensions of the largest stitched image
  unsigned int numOfTiles = m_SettingValue[0] * m_SettingValue[1] * m_SettingValue[2];

  double tileNumber[3];
  tileNumber[0] = m_SettingValue[0];
  tileNumber[1] = m_SettingValue[1];
  tileNumber[2] = m_SettingValue[2];

  double tileSize[3];
  tileSize[0] = m_SettingValue[9];
  tileSize[1] = m_SettingValue[10];
  tileSize[2] = m_SettingValue[11];

  // Read second two lines
  StringVectorType m_TileInfoName;
  m_TileInfoName.reserve( 100 );

  std::getline ( infile, line);
  std::stringstream tileInfoNameStream( line );

  for( unsigned int i = 0; i < 9; i++ )
  {
    std::getline ( tileInfoNameStream, value, ',' );
    m_TileInfoName[i] = value;
    //std::cout << value << std::endl;
  }

  vnlMatrixType m_TileInfoValue (numOfTiles, 12);
  for( unsigned int i = 0; i < numOfTiles; )
  {
    //std::cout << "i = " << i;
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
      for( unsigned int j = 6, k = 0; j < 9; j++, k++ )
      {
        m_TileInfoValue[i][j+3] = m_TileInfoValue[i][j] + tileSize[k];
      }
      //std::cout << std::endl;
      i++;
    }
    else
    {
      //std::cout << std::endl;
    }
  }

  double minStart[3];
  minStart[0] = m_TileInfoValue.get_column(6).min_value();
  minStart[1] = m_TileInfoValue.get_column(7).min_value();
  minStart[2] = m_TileInfoValue.get_column(8).min_value();

  double maxEnd[3];
  maxEnd[0] = m_TileInfoValue.get_column(9).max_value();
  maxEnd[1] = m_TileInfoValue.get_column(10).max_value();
  maxEnd[2] = m_TileInfoValue.get_column(11).max_value();

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
        searchStringXYZT << std::setfill( '0' ) << std::setw( 4 ) << argv[4] << "t.tif";

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
        }
        searchStringCH.str( std::string() );
        searchStringXYZT.str( std::string() );
      }
    }
  }

  // Read one image to get tilePixelDimensions and spacing
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( tileFileNameArray[0][0][0] );
  reader->Update();

  ImageType::Pointer currentImage = reader->GetOutput();
  ImageType::SizeType tilePixelDimension = currentImage->GetLargestPossibleRegion().GetSize();

  ImageType::SpacingType spacing;
  spacing[0] = tileSize[0]/tilePixelDimension[1];
  spacing[1] = tileSize[1]/tilePixelDimension[0];
  spacing[2] = tileSize[2]/tilePixelDimension[2];

  std::cout << tileNumber[0] << ' ' << tileNumber[1] << ' ' << tileNumber[2] << std::endl;
  std::cout << minStart[0] << ' ' << minStart[1] << ' ' << minStart[2] << std::endl;
  std::cout << maxEnd[0] << ' ' << maxEnd[1] << ' ' << maxEnd[2] << std::endl;

  std::cout << tileSize[0] << ' ' << tileSize[1] << ' ' << tileSize[2] << std::endl;
  std::cout << tilePixelDimension << std::endl;
  std::cout << spacing << std::endl;



  // Create the dimensions of the large image


  return EXIT_SUCCESS;
  }
