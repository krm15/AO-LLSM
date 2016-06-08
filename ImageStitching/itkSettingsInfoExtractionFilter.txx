/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 1658 $  // Revision of last commit
  Date: $Date: 2010-06-14 15:49:25 -0400 (Mon, 14 Jun 2010) $  // Date of last commit
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

#ifndef __itkSettingsInfoExtractionFilter_txx
#define __itkSettingsInfoExtractionFilter_txx

#include "itkSettingsInfoExtractionFilter.h"

namespace itk
{
template < class TSegmentImage >
SettingsInfoExtractionFilter< TSegmentImage >
::SettingsInfoExtractionFilter()
{
  m_Dimension = 3;
  m_SettingName.resize( 100 );
  m_SettingValue.resize( 100 );

  m_NumberOfTiles = 1;
}


template < class TSegmentImage >
void
SettingsInfoExtractionFilter< TSegmentImage >::
Read( std::istream& os )
{
  // Read first two lines
  std::string value, line;
  std::getline ( os, line);
  std::stringstream nameStream(line);

  std::getline ( os, line);
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

  for( unsigned int i = 0; i < m_Dimension; i++ )
  {
    m_NumberOfTiles *= m_SettingValue[i];
    m_TileNumber[i] = m_SettingValue[i];
    m_TileSize[i] = m_SettingValue[9+i];
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

}

} /* end namespace itk */

#endif
