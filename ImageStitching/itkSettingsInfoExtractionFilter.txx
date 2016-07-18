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

static itk::SimpleFastMutexLock m_MutexSIEF;

namespace itk
{
template < class TValueType, class TInputImage >
SettingsInfoExtractionFilter< TValueType, TInputImage >
::SettingsInfoExtractionFilter()
{
  m_SettingFieldName.resize( 100 );
  m_SettingFieldValue.resize( 100 );

  m_NumberOfTiles = 1;
  m_StitchedImage = ITK_NULLPTR;
  m_OffsetFile = "";
  m_ChannelPrefix = "_ch";
  m_RegisterZTiles = false;
  m_SharedData = ITK_NULLPTR;
  m_ZTile = 0;
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
UpdateTileCoverage( std::istream& os )
{
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    for( unsigned int j = 0; j < 2; j++ )
    {
      for( unsigned int k = 0; k < 2; k++ )
      {
        m_SharedData->m_TileCover[i][j][k].resize( m_TileNumber[i] );
      }
    }
  }

  for( unsigned int j = 0; j < ImageDimension; j++ )
  {
    for( unsigned int i = 0; i < m_NumberOfTiles; i++ )
    {
      unsigned int temp =  m_TileInfoValue[i][j];
      m_SharedData->m_TileCover[j][0][0][temp] = m_TransformedTileInfoValue[i][j+3];
      m_SharedData->m_TileCover[j][1][0][temp] = m_TransformedTileInfoValue[i][j+6];
    }
  }

  for( unsigned int j = 0; j < ImageDimension; j++ )
  {
    m_MinimumStart[j] = 1000000.0;
    m_MaximumEnd[j]   = -m_MinimumStart[j];
    for( unsigned int k = 0; k < m_TileNumber[j]; k++ )
    {
      // Clipping should be based on offsets as well
      m_SharedData->m_TileCover[j][0][1][k] = m_SharedData->m_TileCover[j][0][0][k]
          + 0.5*m_TileOverlap[j];
      m_SharedData->m_TileCover[j][1][1][k] = m_SharedData->m_TileCover[j][1][0][k]
          - 0.5*m_TileOverlap[j];

      if ( m_MinimumStart[j] > m_SharedData->m_TileCover[j][0][0][k] )
      {
        m_MinimumStart[j] = m_SharedData->m_TileCover[j][0][0][k];
      }

      if ( m_MaximumEnd[j] < m_SharedData->m_TileCover[j][1][0][k] )
      {
        m_MaximumEnd[j] = m_SharedData->m_TileCover[j][1][0][k];
      }
    }
  }

  // Identify first and last tile in each dimension
  // Extend the coverage of the first and last tile by overlap
  for( unsigned int j = 0; j < ImageDimension; j++ )
  {
    unsigned int firstTile(0), lastTile(m_TileNumber[j]-1);

    for( unsigned int k = 0; k < m_TileNumber[j]; k++ )
    {
      if ( m_SharedData->m_TileCover[j][0][1][firstTile] > m_SharedData->m_TileCover[j][0][1][k] )
      {
        firstTile = k;
      }
      if ( m_SharedData->m_TileCover[j][1][1][lastTile] < m_SharedData->m_TileCover[j][1][1][k] )
      {
        lastTile = k;
      }
    }
    m_SharedData->m_TileCover[j][0][1][firstTile] = m_MinimumStart[j];
    m_SharedData->m_TileCover[j][1][1][lastTile] = m_MaximumEnd[j];
  }
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
RegisterTiles()
{
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_TileOffset[i].resize( m_TileNumber[2], 0.0 );
    m_TileEffectiveOffset[i].resize( m_TileNumber[2], 0.0 );
  }

  // Read offset file
  if ( !m_OffsetFile.empty() )
  {
    ReadOffsetFile();
  }

  if ( !m_RegisterZTiles )
  {
    return;
  }

  // Step through z tiles
  unsigned int i = m_ZTile;

  PointType origin;
  origin[0] = m_MinimumStart[0];
  origin[1] = m_MinimumStart[1];
  origin[2] = m_SharedData->m_TileCover[2][1][0][m_ZTile] - m_TileOverlap[2];

  IndexType  roiIndex;
  roiIndex.Fill( 0 );

  SizeType roiSize = m_StitchDimension;
  roiSize[2] = static_cast<SizeValueType>( m_TileOverlap[2]/m_TileSpacing[2] ) + 1;

  RegionType roi;
  roi.SetIndex( roiIndex );
  roi.SetSize( roiSize );

  // Assemble ROI of static and moving images
  ImagePointer m_ROIStaticImage = ImageType::New();
  m_ROIStaticImage->SetOrigin( origin );
  m_ROIStaticImage->SetSpacing( m_TileSpacing );
  m_ROIStaticImage->SetRegions( roi );
  m_ROIStaticImage->Allocate();
  m_ROIStaticImage->FillBuffer( 0.0 );

  FillROIFilterPointer fillStaticROI = FillROIFilterType::New();
  fillStaticROI->SetInput( m_ROIStaticImage );
  fillStaticROI->SetZTile( m_ZTile );
  fillStaticROI->SetSharedData( m_SharedData );
  fillStaticROI->InPlaceOn();
  fillStaticROI->SetNumberOfThreads( 1 );
  fillStaticROI->Update();
  ImagePointer m_staticImage = fillStaticROI->GetOutput();
  std::cout << "Extracted static image" << std::endl;

  WriterPointer writer1 = WriterType::New();
  writer1->SetInput( m_staticImage );
  writer1->SetFileName( "Static.mha" );
  writer1->Update();

  ImagePointer m_ROIMovingImage = ImageType::New();
  m_ROIMovingImage->SetOrigin( origin );
  m_ROIMovingImage->SetSpacing( m_TileSpacing );
  m_ROIMovingImage->SetRegions( roi );
  m_ROIMovingImage->Allocate();
  m_ROIMovingImage->FillBuffer( 0.0 );

  FillROIFilterPointer fillMovingROI = FillROIFilterType::New();
  fillMovingROI->SetInput( m_ROIMovingImage );
  fillMovingROI->SetZTile( m_ZTile );
  fillMovingROI->SetSharedData( m_SharedData );
  fillMovingROI->InPlaceOn();
  fillMovingROI->SetNumberOfThreads( 1 );
  fillMovingROI->Update();
  ImagePointer m_movingImage = fillMovingROI->GetOutput();
  std::cout << "Extracted moving image" << std::endl;

  WriterPointer writer2 = WriterType::New();
  writer2->SetInput( m_movingImage );
  writer2->SetFileName( "Moving.mha" );
  writer2->Update();

  // Compute overlap of ROI

  MetricPointer         metric        = MetricType::New();
  TransformPointer      transform     = TransformType::New();
  OptimizerPointer      optimizer     = OptimizerType::New();
  optimizer->SetMaximumStepLength( 3.00 );
  optimizer->SetMinimumStepLength( 0.01 );
  optimizer->SetNumberOfIterations( 1 );
  InterpolatorPointer   interpolator  = InterpolatorType::New();
  ParametersType initialParameters( transform->GetNumberOfParameters() );

  double bestValue = 1000000, besti, bestj, bestk, value;
  value = bestValue;
  for( float i = -2.0; i <= 2.0; i+=0.1 )
  {
    for( float j = -2.0; j <= 2.0; j+=0.1 )
    {
      for( float k = -2.0; k <= 2.0; k+=0.2 )
      {
        std::cout << i<< ' ' << j << ' ' << k << ' ' << value << std::endl;
        initialParameters[0] = m_TileOffset[0][m_ZTile] + i;  // Initial offset along X
        initialParameters[1] = m_TileOffset[1][m_ZTile] + j;  // Initial offset along Y
        initialParameters[2] = m_TileOffset[2][m_ZTile] + k;  // Initial offset along Z

        RegistrationPointer registration  = RegistrationType::New();
        registration->SetNumberOfThreads( 1 );
        registration->SetMetric(        metric        );
        registration->SetOptimizer(     optimizer     );
        registration->SetTransform(     transform     );
        registration->SetInterpolator(  interpolator  );
        registration->SetFixedImage(    m_staticImage    );
        registration->SetMovingImage(   m_movingImage   );
        registration->SetFixedImageRegion( m_staticImage->GetLargestPossibleRegion() );
        registration->SetInitialTransformParameters( initialParameters );
        registration->Update();

        value = static_cast<double>( optimizer->GetValue() );
        if ( value  < bestValue )
        {
          bestValue = value;
          besti = m_TileOffset[0][m_ZTile] + i;
          bestj = m_TileOffset[1][m_ZTile] + j;
          bestk = m_TileOffset[2][m_ZTile] + k;
          std::cout << besti<< ' ' << bestj << ' ' << bestk << ' ' << bestValue << std::endl;
        }
      }
    }
  }
  std::cout << besti<< ' ' << bestj << ' ' << bestk << ' ' << bestValue << std::endl;
  m_TileOffset[0][m_ZTile] = besti;
  m_TileOffset[1][m_ZTile] = bestj;
  m_TileOffset[2][m_ZTile] = bestk;

  // Write out the offsets
  if ( !m_OffsetFile.empty() )
  {
    WriteOffsetFile();
  }
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
ReadOffsetFile()
{
  std::ifstream os ( m_OffsetFile.c_str() );

  if ( !os )
  {
    return;
  }

  std::string line, value;
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    std::getline ( os, line );
    std::stringstream valueStream( line );
    for( unsigned int j = 0; j < m_TileNumber[2]; j++ )
    {
      std::getline ( valueStream, value, ' ' );
      m_TileOffset[i][j] = atof( value.c_str() );
    }
  }

  // Compile offsets to identify effective offsets
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_TileEffectiveOffset[i][0] = 0.0;
    for( unsigned int j = 1; j < m_TileNumber[2]; j++ )
    {
      m_TileEffectiveOffset[i][j] = m_TileEffectiveOffset[i][j-1] + m_TileOffset[i][j];
    }
  }
  os.close();
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
WriteOffsetFile()
{
  std::ofstream os ( m_OffsetFile.c_str() );

  if ( !os )
  {
    std::cout << "error in offset file opening" << std::endl;
    return;
  }

  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    for( unsigned int j = 0; j < m_TileNumber[2]; j++ )
    {
      os << m_TileOffset[i][j] << ' ';
    }
    os << std::endl;
  }

  os.close();
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
ReadTileInfo( std::istream& os )
{
  std::string value, line;

  // Do a string search to determine the locations
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_NumberOfTiles *= m_SettingFieldValue[i];
    m_TileNumber[i] = m_SettingFieldValue[i];
    m_TileSize[i] = m_SettingFieldValue[9+i];
    m_TileOverlap[i] = m_SettingFieldValue[6+i];
  }

  // Read next two lines
  StringVectorType m_TileInfoName;
  m_TileInfoName.resize( 100 );

  std::getline ( os, line);
  std::stringstream tileInfoNameStream( line );

  for( unsigned int i = 0; i < 9; i++ )
  {
    std::getline ( tileInfoNameStream, value, ',' );
    m_TileInfoName[i] = value;
    //std::cout << value << std::endl;
  }

  m_TileInfoValue.set_size( m_NumberOfTiles, 9 );
  for( unsigned int i = 0; i < m_NumberOfTiles; )
  {
    //    std::cout << i << std::endl;
    std::getline ( os, line);

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
  }
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
TransformCoordinateAxes()
{
  vnlVectorType tileCenter, stitchCenter, newCenter;
  tileCenter.set_size( ImageDimension );
  stitchCenter.set_size( ImageDimension );
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    stitchCenter[i] = m_TileInfoValue.get_column(i+6).mean();
    //std::cout << stitchCenter[i] << ' ';
  }

  double theta = 31.8 * vnl_math::pi_over_180;
  vnlMatrixType rotMatrix( ImageDimension, ImageDimension );
  rotMatrix[0][0] = rotMatrix[0][2] = rotMatrix[1][1] = rotMatrix[2][1] = 0.0;
  rotMatrix[0][1] = -1;
  rotMatrix[1][0] = -vcl_cos( theta );
  rotMatrix[2][2] = vcl_cos( theta );
  rotMatrix[1][2] = rotMatrix[2][0] = vcl_sin( theta );

  //std::cout << rotMatrix << std::endl;

  m_TransformedTileInfoValue.set_size( m_NumberOfTiles, 3*ImageDimension );
  for( unsigned int i = 0; i < m_NumberOfTiles; i++ )
  {
    for( unsigned int j = 0; j < ImageDimension; j++ )
    {
      tileCenter[j] = m_TileInfoValue[i][j+6];
      //std::cout << m_TileInfoValue[i][j] << ' ';
    }

    newCenter = rotMatrix * ( tileCenter - stitchCenter );

    //std::cout << newCenter << std::endl;

    for( unsigned int j = 0; j < ImageDimension; j++ )
    {
      m_TransformedTileInfoValue[i][j] = newCenter[j];
      m_TransformedTileInfoValue[i][j+ImageDimension] = newCenter[j] - 0.5*m_TileSize[j];
      m_TransformedTileInfoValue[i][j+2*ImageDimension] = newCenter[j] + 0.5*m_TileSize[j];
    }
  }
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
UpdateFileNameLookup( std::istream& os )
{
  // Read all the files in the input directory of type ch and at timePoint
  std::string filename;
  std::stringstream searchStringCH, searchStringXYZT, filename2;

  //Identify the channel type
  bool ChannelNameSet = false;
  std::string searchString = "nm";

  DirectoryPointer directory = DirectoryType::New();
  directory->Load( m_TileDirectory.c_str() );

  m_SharedData->m_TileFileNameArray.resize( m_SettingFieldValue[0] );
  unsigned int m_TrueCountOfTiles = 0;
  for( unsigned int i = 0; i < m_TileNumber[0]; i++ )
  {
    m_SharedData->m_TileFileNameArray[i].resize( m_TileNumber[1] );
    for( unsigned int j = 0; j < m_TileNumber[1]; j++ )
    {
      m_SharedData->m_TileFileNameArray[i][j].resize( m_TileNumber[2] );
      for( unsigned int k = 0; k < m_TileNumber[2]; k++ )
      {
        m_SharedData->m_TileFileNameArray[i][j][k] = std::string();
        searchStringCH << m_ChannelPrefix << m_ChannelNumber;
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << i << "x_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << j << "y_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << k << "z_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 4 ) << m_TimePoint << "t";

        for ( unsigned int m = 0; m < directory->GetNumberOfFiles(); m++)
        {
          filename = directory->GetFile( m );

          if ( ( filename.find( searchStringCH.str() ) != std::string::npos ) &&
               ( filename.find( searchStringXYZT.str() ) != std::string::npos ) )
          {
            //std::cout << i << ' ' << j << ' ' << k << ' ' << filename << std::endl;
            filename2 << m_TileDirectory << filename;
            m_SharedData->m_TileFileNameArray[i][j][k] = filename2.str();

            if ( !ChannelNameSet )
            {
              unsigned int pos = filename.find( searchString );
              m_ChannelName = filename.substr( pos-3, 5 );
              ChannelNameSet = true;
              m_SampleName = filename2.str();
            }
            m_TrueCountOfTiles++;
          }
          filename2.str( std::string() );
        }
        searchStringCH.str( std::string() );
        searchStringXYZT.str( std::string() );
      }
    }
  }
  m_NumberOfTiles = m_TrueCountOfTiles;
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
Read()
{
  if ( !m_SharedData )
  {
    m_SharedData = SharedDataType::New();
  }

  DirectoryPointer directory = DirectoryType::New();
  directory->Load( m_SettingsDirectory.c_str() );

  bool foundFile = false;
  std::string filename;
  for ( unsigned m = 0; m < directory->GetNumberOfFiles(); m++)
  {
    filename = directory->GetFile( m );

    if ( ( filename.find( "3D settings" ) != std::string::npos ) &&
         ( filename.find( ".csv" ) != std::string::npos ) &&
         ( !foundFile ) )
    {
      foundFile = true;
      break;
    }
  }

  if ( !foundFile )
  {
    std::cout << "3D settings file not found" << std::endl;
    return;
  }

  filename =  m_SettingsDirectory + filename;

  std::ifstream os ( filename.c_str() );

  if ( !os )
  {
    std::cout << "error in file opening" << std::endl;
    return;
  }

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
    m_SettingFieldName[i] = value;

    //std::cout << value << ' ';
    std::getline ( valueStream, value, ',' );
    m_SettingFieldValue[i] = atof( value.c_str() );
    //std::cout << value << std::endl;
  }

  ReadTileInfo( os );
  std::cout << "Read tile info" << std::endl;

  // Create a lookup of filenames
  UpdateFileNameLookup( os );
  std::cout << "Updated file name lookup" << std::endl;

  TransformCoordinateAxes();
  std::cout << "Transformed coordinate axes" << std::endl;

  // Create a vector of tile origins along each axis for given timepoint
  UpdateTileCoverage( os );
  std::cout << "Updated tile coverage" << std::endl;

  // Read one image to get m_TileDimensions and m_TileSpacing
  {
    ReaderPointer reader = ReaderType::New();
    reader->SetFileName ( m_SampleName );
    reader->Update();
    ImagePointer currentImage = reader->GetOutput();
    m_TileDimension = currentImage->GetLargestPossibleRegion().GetSize();

    m_TileSpacing[0] = m_TileSize[1]/m_TileDimension[0];
    m_TileSpacing[1] = m_TileSize[0]/m_TileDimension[1];
    m_TileSpacing[2] = m_TileSize[2]/m_TileDimension[2];
    std::cout << "Read tile dimensions" << std::endl;
  }

  // Store m_TileDimension correctly since image is flipped
  unsigned int temp = m_TileDimension[0];
  m_TileDimension[0] = m_TileDimension[1];
  m_TileDimension[1] = temp;

  // Create a stitched image
  CreateStitchedImage();

  // Register files
  RegisterTiles();

  os.close();
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
CreateStitchedImage()
{
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_StitchIndex[i]     = 0;
    m_StitchSize[i]      = m_MaximumEnd[i] - m_MinimumStart[i];
    m_StitchOrigin[i]    = m_MinimumStart[i];
    m_StitchDimension[i] = m_StitchSize[i]/m_TileSpacing[i];
  }

  m_StitchedImage = ImageType::New();
  m_StitchedImage->SetOrigin( m_StitchOrigin );
  m_StitchedImage->SetSpacing( m_TileSpacing );
  m_StitchedImage->SetRegions( m_StitchRegion );
}

} /* end namespace itk */

#endif
