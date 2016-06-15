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
template < class TValueType, class TInputImage >
SettingsInfoExtractionFilter< TValueType, TInputImage >
::SettingsInfoExtractionFilter()
{
  m_SettingFieldName.resize( 100 );
  m_SettingFieldValue.resize( 100 );

  m_Blending = false;
  m_NumberOfTiles = 1;
  m_StitchedImage = ITK_NULLPTR;
  m_CorrectionImage = ITK_NULLPTR;
  m_CorrectionDirectory = "";
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
UpdateTileCoverage( std::istream& os )
{
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_TileCoverStart[i].resize( m_TileNumber[i] );
    m_TileCoverEnd[i].resize( m_TileNumber[i] );
  }

  for( unsigned int j = 0; j < ImageDimension; j++ )
  {
    for( unsigned int i = 0; i < m_NumberOfTiles; i++ )
    {
      unsigned int temp =  m_TileInfoValue[i][j];
      m_TileCoverStart[j][temp] = m_TransformedTileInfoValue[i][j+3] + 0.5*m_TileOverlap[j];
      m_TileCoverEnd[j][temp] = m_TransformedTileInfoValue[i][j+6] - 0.5*m_TileOverlap[j];
    }
  }

  // Identify first and last tile in each dimension
  // Extend the coverage of the first and last tile by overlap
  for( unsigned int j = 0; j < ImageDimension; j++ )
  {
    unsigned int firstTile, lastTile;
    firstTile = 0;
    lastTile = m_TileNumber[j]-1;

    double start, stop;
    for( unsigned int k = 0; k < m_TileNumber[j]; k++ )
    {
      if ( m_TileCoverStart[j][firstTile] > m_TileCoverStart[j][k] )
      {
        firstTile = k;
      }
      if ( m_TileCoverEnd[j][lastTile] < m_TileCoverEnd[j][k] )
      {
        lastTile = k;
      }
    }
    m_TileCoverStart[j][firstTile] = m_MinimumStart[j];
    m_TileCoverEnd[j][lastTile] = m_MaximumEnd[j];
  }
}

template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
UpdateStitchDimensions( std::istream& os )
{
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_MinimumStart[i] = m_TransformedTileInfoValue.get_column(i+3).min_value();
    m_MaximumEnd[i] = m_TransformedTileInfoValue.get_column(i+6).max_value();
  }
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
ReadTileInfo( std::istream& os )
{
  std::string value, line;

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
    else
    {
      //std::cout << std::endl;
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

  m_TransformedTileInfoValue.set_size(m_NumberOfTiles, 3*ImageDimension);
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

  m_TileFileNameArray.resize( m_SettingFieldValue[0] );

  for( unsigned int i = 0; i < m_TileNumber[0]; i++ )
  {
    m_TileFileNameArray[i].resize( m_TileNumber[1] );
    for( unsigned int j = 0; j < m_TileNumber[1]; j++ )
    {
      m_TileFileNameArray[i][j].resize( m_TileNumber[2] );
      for( unsigned int k = 0; k < m_TileNumber[2]; k++ )
      {
        m_TileFileNameArray[i][j][k] = std::string();
        searchStringCH << "_ch" << m_ChannelNumber;
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << i << "x_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << j << "y_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << k << "z_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 4 ) << m_TimePoint << "t";

        //std::cout << i << j << k << std::endl;
        for ( unsigned m = 0; m < directory->GetNumberOfFiles(); m++)
        {
          filename = directory->GetFile( m );

          if ( ( filename.find( searchStringCH.str() ) != std::string::npos ) &&
               ( filename.find( searchStringXYZT.str() ) != std::string::npos ) )
          {
            //std::cout << filename << std::endl;
            filename2 << m_TileDirectory << filename;
            m_TileFileNameArray[i][j][k] = filename2.str();

            if ( !ChannelNameSet )
            {
              unsigned int pos = filename.find( searchString );
              m_ChannelName = filename.substr( pos-3, 5 );
              ChannelNameSet = true;
              m_SampleName = filename2.str();
            }

          }
          filename2.str( std::string() );
        }
        searchStringCH.str( std::string() );
        searchStringXYZT.str( std::string() );
      }
    }
  }
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
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
    m_SettingFieldName[i] = value;

    //std::cout << value << ' ';
    std::getline ( valueStream, value, ',' );
    m_SettingFieldValue[i] = atof( value.c_str() );
    //std::cout << value << std::endl;
  }

  ReadTileInfo( os );
  std::cout << "Read tile info" << std::endl;

  TransformCoordinateAxes();
  std::cout << "Transformed coordinate axes" << std::endl;

  // Identify total tile coverage
  UpdateStitchDimensions( os );
  std::cout << "Updated stitch dimensions" << std::endl;

  // Create a vector of tile origins along each axis for given timepoint
  UpdateTileCoverage( os );
  std::cout << "Updated tile coverage" << std::endl;

  // Create a lookup of filenames
  UpdateFileNameLookup( os );
  std::cout << "Updated file name lookup" << std::endl;

  // Read one image to get m_TileDimensions and m_TileSpacing
  ReaderPointer reader = ReaderType::New();
  reader->SetFileName ( m_SampleName );
  reader->SetGlobalWarningDisplay( 0 );
  reader->Update();
  ImagePointer currentImage = reader->GetOutput();
  m_TileDimension = currentImage->GetLargestPossibleRegion().GetSize();

  m_TileSpacing[0] = m_TileSize[1]/m_TileDimension[0];
  m_TileSpacing[1] = m_TileSize[0]/m_TileDimension[1];
  m_TileSpacing[2] = m_TileSize[2]/m_TileDimension[2];
  std::cout << "Read tile dimensions" << std::endl;

  // Read the correction image
  ReadCorrectionImage();
  std::cout << "Read correction image" << std::endl;
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
ReadCorrectionImage()
{
  if ( m_CorrectionDirectory.empty() )
  {
    return;
  }

  bool IsCorrectionImage = false;
  std::string filename;
  std::stringstream filename2;
  double variance = 1.0;
  double threshold = 30.0;

  RImagePointer currentImage;

  DirectoryPointer directory = DirectoryType::New();
  directory->Load( m_CorrectionDirectory.c_str() );
  for ( unsigned int m = 0; m < directory->GetNumberOfFiles(); m++)
  {
    filename = directory->GetFile( m );
    if ( ( ! filename.empty() ) && ( !IsCorrectionImage ) )
    {
      if ( filename.find( m_ChannelName ) != std::string::npos )
      {
        filename2 << m_CorrectionDirectory.c_str() << filename;

        // Read the correction image
        RReaderPointer reader = RReaderType::New();
        reader->SetFileName ( filename2.str() );
        reader->SetGlobalWarningDisplay( 0 );
        reader->Update();

        GaussianFilterPointer gaussianFilter = GaussianFilterType::New();
        gaussianFilter->SetInput( reader->GetOutput() );
        gaussianFilter->SetVariance( variance );
        gaussianFilter->Update();

        currentImage = gaussianFilter->GetOutput();
        currentImage->DisconnectPipeline();
        break;
      }
      filename2.str( std::string() );
    }
  }

  RSpacingType sp;
  sp[0] = m_TileSpacing[0];
  sp[1] = m_TileSpacing[1];
  currentImage->SetSpacing( sp );

  RRegionType rroi;

  RSizeType rsize = currentImage->GetLargestPossibleRegion().GetSize();

  RIndexType rindex;
  rindex[0] = 0.5*( rsize[0] - m_TileDimension[0] );
  rindex[1] = 0.5*( rsize[1] - m_TileDimension[0] );
  rroi.SetIndex( rindex );

  rsize[0] = m_TileDimension[0];
  rsize[1] = m_TileDimension[1];
  rroi.SetSize( rsize );

  ROIFilterPointer roiFilter = ROIFilterType::New();
  roiFilter->SetRegionOfInterest( rroi );
  roiFilter->SetInput( currentImage );
  roiFilter->Update();
  m_CorrectionImage = roiFilter->GetOutput();
  m_CorrectionImage->DisconnectPipeline();

  RIteratorType It( m_CorrectionImage, m_CorrectionImage->GetLargestPossibleRegion() );
  It.GoToBegin();
  while( !It.IsAtEnd() )
  {
    double p = It.Get();
    if ( p < threshold )
    {
      It.Set( threshold );
    }
    ++It;
  }
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


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
AllocateROI()
{
  if ( !m_StitchedImage )
  {
    std::cout <<  "Stitch image created" << std::endl;
    CreateStitchedImage();
  }

  m_ROIImage = ImageType::New();
  m_ROIImage->SetOrigin( m_ROIOrigin );
  m_ROIImage->SetSpacing( m_TileSpacing );
  m_ROIImage->SetRegions( m_ROI );
  m_ROIImage->Allocate();
  m_ROIImage->FillBuffer( 0.0 );

  // Identify all the tiles that belong to this roi
  std::cout << "Setting scan start and end values for ROI" << std::endl;
  for( unsigned int k = 0; k < ImageDimension; k++ )
  {
    double beginCorner = m_ROIOrigin[k];
    double endCorner = m_ROIOrigin[k] + m_ROI.GetSize()[k] * m_TileSpacing[k];

    //std::cout <<  beginCorner << ' ' << endCorner << std::endl;
    double scanStartVal = 100000, scanEndVal = -100000;
    for( unsigned int i = 0; i < m_TileNumber[k]; i++ )
    {
      //std::cout << m_TileCoverStart[k][i] << ' '
      //          << m_TileCoverEnd[k][i] << std::endl;
      if ( ( beginCorner >= m_TileCoverStart[k][i] - 0.0001 ) &&
           ( beginCorner <= m_TileCoverEnd[k][i] + 0.0001 ) &&
           ( scanStartVal >= m_TileCoverStart[k][i] ) )
      {
        m_ScanStart[k] = i;
        scanStartVal =  m_TileCoverStart[k][i];
      }

      if ( ( endCorner >= m_TileCoverStart[k][i] - 0.001 ) &&
           ( endCorner <= m_TileCoverEnd[k][i] + 0.001 ) &&
           ( scanEndVal <= m_TileCoverEnd[k][i] ) )
      {
        m_ScanEnd[k] = i;
        scanEndVal =  m_TileCoverEnd[k][i];
      }
    }

    unsigned int temp;
    if ( m_ScanStart[k] > m_ScanEnd[k] )
    {
      temp = m_ScanEnd[k];
      m_ScanEnd[k] = m_ScanStart[k];
      m_ScanStart[k] = temp;
    }

    //std::cout << m_ScanStart[k] << ' ' << m_ScanEnd[k] << std::endl;
  }

  FillROI();
  std::cout << "ROI filled" << std::endl;
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >
::OverlapRegion( ImagePointer A, ImagePointer B,
  RegionType& rA, RegionType& rB )
{
  SizeType sizeA, sizeB, s;
  sizeA = A->GetLargestPossibleRegion().GetSize();
  sizeB = B->GetLargestPossibleRegion().GetSize();

  IndexType sIndexA, sIndexB;
  IndexType tIndexA, tIndexB;

  A->TransformPhysicalPointToIndex( B->GetOrigin(), tIndexA );
  B->TransformPhysicalPointToIndex( A->GetOrigin(), tIndexB );

  PointType originA = A->GetOrigin();
  PointType originB = B->GetOrigin();

  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    if ( originA[i] > originB[i] )
    {
      sIndexA[i] = 0;
      sIndexB[i] = tIndexB[i];
      s[i] = sizeA[i];
      if ( s[i] > static_cast< SizeValueType >( sizeB[i] - sIndexB[i] - 1 ) )
      {
        s[i] = sizeB[i] - sIndexB[i];
      }
    }
    else
    {
      sIndexB[i] = 0;
      sIndexA[i] = tIndexA[i];
      s[i] = sizeB[i];
      if ( s[i] > static_cast< SizeValueType >(
        sizeA[i] - sIndexA[i] - 1 ) )
      {
        s[i] = sizeA[i] - sIndexA[i];
      }
    }
  }

  rA.SetIndex( sIndexA );
  rA.SetSize( s );
  rB.SetIndex( sIndexB );
  rB.SetSize( s );
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
FillROI()
{
  FixedArray<unsigned int, 3> axesOrder;
  axesOrder[0] = 1;
  axesOrder[1] = 0;
  axesOrder[2] = 2;

  PixelType p;

  // Start a loop that will read all the tiles from zScanStart to zScanEnd
  PointType currentTileOrigin;
  RegionType currentTileRegion, roiSubRegion;
  IndexType temp;
  for( unsigned int i = m_ScanStart[0]; i <= m_ScanEnd[0]; i++ )
  {
    currentTileOrigin[0] = m_TileCoverStart[0][i];
    for( unsigned int j = m_ScanStart[1]; j <= m_ScanEnd[1]; j++ )
    {
      currentTileOrigin[1] = m_TileCoverStart[1][j];
      for( unsigned int k = m_ScanStart[2]; k <= m_ScanEnd[2]; k++ )
      {
        currentTileOrigin[2] = m_TileCoverStart[2][k];

        m_ROIImage->TransformPhysicalPointToIndex( currentTileOrigin, temp );
        std::string filename = m_TileFileNameArray[i][j][k];

        if ( m_ROI.IsInside( temp ) && ( ! filename.empty() )  )
        {
          //std::cout << i << ' ' << j << ' ' << k << std::endl;
          //std::cout << filename.c_str() << std::endl;

          ReaderPointer reader = ReaderType::New();
          reader->SetFileName( filename.c_str() );
          reader->SetGlobalWarningDisplay( 0 );
          reader->Update();
          ImagePointer cImage = reader->GetOutput();
          cImage->DisconnectPipeline();

          if ( m_CorrectionImage )
          {
            IteratorType cIt( cImage, cImage->GetLargestPossibleRegion() );
            cIt.GoToBegin();
            RIteratorType corrIt( m_CorrectionImage, m_CorrectionImage->GetLargestPossibleRegion() );
            corrIt.GoToBegin();
            while( !cIt.IsAtEnd() )
            {
              if ( corrIt.IsAtEnd() )
              {
                corrIt.GoToBegin();
              }
              p = static_cast<PixelType>( 100*double(cIt.Get())/corrIt.Get() );
              cIt.Set( p );
              ++cIt;
              ++corrIt;
            }
          }

          PermuteAxesFilterPointer pAFilter = PermuteAxesFilterType::New();
          pAFilter->SetInput( cImage );
          pAFilter->SetOrder( axesOrder );
          pAFilter->Update();
          ImagePointer pImage = pAFilter->GetOutput();

          ImagePointer currentImage = ImageType::New();
          currentImage->SetOrigin( pImage->GetOrigin() );
          currentImage->SetSpacing( pImage->GetSpacing() );
          currentImage->SetRegions( pImage->GetLargestPossibleRegion() );
          currentImage->Allocate();

          IteratorType pIt( pImage, pImage->GetLargestPossibleRegion() );
          IteratorType currentIt( currentImage, currentImage->GetLargestPossibleRegion() );
          while( !pIt.IsAtEnd() )
          {
            currentIt.Set( pIt.Get() );
            ++pIt;
            ++currentIt;
          }

          currentImage->SetOrigin( currentTileOrigin );

          //std::cout << "Current tile origin" << std::endl;
          //std::cout << currentTileOrigin << std::endl;

          OverlapRegion( currentImage, m_ROIImage, currentTileRegion, roiSubRegion );

          //std::cout << "Current tile region" << std::endl;
          //std::cout << currentTileRegion << std::endl;

          //std::cout << "ROI region" << std::endl;
          //std::cout << roiSubRegion << std::endl;


          // Using these images, fill up roiImage
          IteratorType rIt( m_ROIImage, roiSubRegion );
          rIt.GoToBegin();

          IteratorType tIt( currentImage, currentTileRegion );
          tIt.GoToBegin();

          while( !tIt.IsAtEnd() )
          {
            rIt.Set( tIt.Get() );
            ++tIt;
            ++rIt;
          }
        }
      }
    }
  }
}

} /* end namespace itk */

#endif
