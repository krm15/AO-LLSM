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
  m_OffsetFilePath = "";
  m_ChannelPrefix = "_ch";
  m_RegisterZTiles = false;
  m_SharedData = ITK_NULLPTR;
  m_ZTileStart = 0;
  m_ZTileEnd = 0;
  //m_StepLength = 0.5;
  //m_SearchRadius = 5.0;
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
OverlapRegion( ImagePointer A, ImagePointer B, RegionType& rA, RegionType& rB )
{
  SizeType sizeA, sizeB, s;
  sizeA = A->GetLargestPossibleRegion().GetSize();
  sizeB = B->GetLargestPossibleRegion().GetSize();

  PointType originA = A->GetOrigin();
  PointType originB = B->GetOrigin();

  IndexType sIndexA, sIndexB;
  IndexType tIndexA, tIndexB;

  A->TransformPhysicalPointToIndex( originB, tIndexA );
  B->TransformPhysicalPointToIndex( originA, tIndexB );

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
      unsigned int temp =  m_TileInfoValue[i][j];// range is m_TileNumber[j]
      m_SharedData->m_TileCover[j][0][0][temp] = m_TransformedTileInfoValue[i][j+3];
      m_SharedData->m_TileCover[j][1][0][temp] = m_TransformedTileInfoValue[i][j+6];
    }
  }

  // HACK
  for( unsigned int j = 0; j < ImageDimension; j++ )
  {
    bool sign = true;
    if ( m_SharedData->m_TileCover[j][0][0][0] > m_SharedData->m_TileCover[j][0][0][1] )
    {
      sign = false;
    }

    for( unsigned int k = 0; k < m_TileNumber[j]; k++ )
    {
       m_SharedData->m_TileCover[j][0][0][k] = ( k * ( m_TileSize[j] -  m_TileOverlap[j]) );

       if (!sign)
       {
         m_SharedData->m_TileCover[j][0][0][k] = -m_SharedData->m_TileCover[j][0][0][k];
       }

       m_SharedData->m_TileCover[j][1][0][k] = m_SharedData->m_TileCover[j][0][0][k] + m_TileSize[j];
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
RegisterTiles( float searchRadius, float stepLength )
{
  bool tileForRegistration = false;
  unsigned int i_,j_;

  if ( m_ZTileEnd > m_TileNumber[2]-2 )
  {
    m_ZTileEnd = m_TileNumber[2]-2;
  }

  for( unsigned int ztile = m_ZTileStart; ztile <= m_ZTileEnd; ztile++ )
  {    
    i_ = m_TileNumber[0]/2;
    j_ = m_TileNumber[1]/2;
    if ( ( !m_SharedData->m_TileFileNameArray[i_][j_][ztile].empty() ) &&
    ( !m_SharedData->m_TileFileNameArray[i_][j_][ztile+1].empty() ) )
    {
      tileForRegistration = true;
    }

    for( unsigned int i = 0; i < m_TileNumber[0]; i++ )
    {
      for( unsigned int j = 0; j < m_TileNumber[1]; j++ )
      {
        if ( ( !m_SharedData->m_TileFileNameArray[i][j][ztile].empty() ) &&
             ( !m_SharedData->m_TileFileNameArray[i][j][ztile+1].empty() ) &&
             ( !tileForRegistration ) )
        {
          tileForRegistration = true;
          i_ = i;
          j_ = j;
        }
      }
    }

    ReaderPointer sreader = ReaderType::New();
    sreader->SetFileName( m_SharedData->m_TileFileNameArray[i_][j_][ztile] );
    sreader->Update();
    ImagePointer staticImage = sreader->GetOutput();
    staticImage->DisconnectPipeline();
    
    ReaderPointer mreader = ReaderType::New();
    mreader->SetFileName( m_SharedData->m_TileFileNameArray[i_][j_][ztile+1] );
    mreader->Update();
    ImagePointer movingImage = mreader->GetOutput();
    movingImage->DisconnectPipeline();

//    std::cout << m_SharedData->m_TileFileNameArray[i_][j_][ztile] << std::endl;
//    std::cout << m_SharedData->m_TileFileNameArray[i_][j_][ztile+1] << std::endl;

    PointType sorigin;
    sorigin[0] = 0.0;
    sorigin[1] = 0.0;
    sorigin[2] = 0.0;

    PointType morigin;
    morigin[0] = m_SharedData->m_TileOffset[0][ztile+1];
    morigin[1] = m_SharedData->m_TileOffset[1][ztile+1];
    morigin[2] = m_TileSize[2] - m_TileOverlap[2]
        + m_SharedData->m_TileOffset[2][ztile+1];

    staticImage->SetOrigin( sorigin );
    movingImage->SetOrigin( morigin );

    staticImage->SetSpacing( m_TileSpacing );
    movingImage->SetSpacing( m_TileSpacing );
    
  //   WriterPointer writer1 = WriterType::New();
  //   writer1->SetInput( staticImage );
  //   writer1->SetFileName( "/home/krm15/output/static.mha" );
  //   writer1->Update();
  // 
  //   WriterPointer writer2 = WriterType::New();
  //   writer2->SetInput( movingImage );
  //   writer2->SetFileName( "/home/krm15/output/moving.mha" );
  //   writer2->Update();  
    
    RegionType sROI, mROI;
    PointType norigin;
    
    double val1(0.0), val2(0.0);
    OverlapRegion( staticImage, movingImage, sROI, mROI );
    IteratorType sIt( staticImage, sROI );
    IteratorType mIt( movingImage, mROI );
    sIt.GoToBegin();
    mIt.GoToBegin();
    while( !sIt.IsAtEnd() )
    {
      val1 += sIt.Get();
      val2 += mIt.Get();
      ++sIt;
      ++mIt;
    }
    double scaleFactor = val1/val2;
    
    double bestValue = std::numeric_limits<double>::max();
    double besti, bestj, bestk, value;
    value = bestValue;

    for( float i = -searchRadius; i <= searchRadius; i+=stepLength )
    {
      norigin[0] = morigin[0] + i;
      for( float j = -searchRadius; j <= searchRadius; j+=stepLength )
      {
        norigin[1] = morigin[1] + j;
        for( float k = -searchRadius; k <= searchRadius; k+=stepLength )
        {
          //std::cout << i<< ' ' << j << ' ' << k << ' ' << value << std::endl;
          norigin[2] = morigin[2] + k;
          movingImage->SetOrigin( norigin );

          OverlapRegion( staticImage, movingImage, sROI, mROI );

          value = 0.0;
          IteratorType sIt( staticImage, sROI );
          IteratorType mIt( movingImage, mROI );
          sIt.GoToBegin();
          mIt.GoToBegin();
          while( !sIt.IsAtEnd() )
          {
            val1 = ( sIt.Get() - scaleFactor * mIt.Get() );
            value += val1 * val1;
            ++sIt;
            ++mIt;
          }
          value /= sROI.GetNumberOfPixels();

          if ( value  < bestValue )
          {
            bestValue = value;
            besti = i;
            bestj = j;
            bestk = k;
            std::cout << "*****" << besti<< ' ' << bestj << ' '
                      << bestk << ' ' << bestValue << std::endl;
          }
        }
      }
    }
    std::cout << besti<< ' ' << bestj << ' ' <<
                 bestk << ' ' << bestValue << std::endl;
    m_SharedData->m_TileOffset[0][ztile+1] += besti;
    m_SharedData->m_TileOffset[1][ztile+1] += bestj;
    m_SharedData->m_TileOffset[2][ztile+1] += bestk;
  }
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
ReadOffsetFile()
{
  for( unsigned int id = 1; id < m_TileNumber[2]; id++  )
  {
    std::stringstream m_OffsetFile;
    m_OffsetFile << m_OffsetFilePath << id << ".txt";
    std::ifstream os ( m_OffsetFile.str().c_str() );

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
        m_SharedData->m_TileOffset[i][j] += atof( value.c_str() );
      }
    }
    os.close();
  }
    
  // Compile offsets to identify effective offsets
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_SharedData->m_TileEffectiveOffset[i][0] = 0.0;
    for( unsigned int j = 1; j < m_TileNumber[2]; j++ )
    {
      m_SharedData->m_TileEffectiveOffset[i][j] =
          m_SharedData->m_TileEffectiveOffset[i][j-1] +
          m_SharedData->m_TileOffset[i][j];
    }
  }
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
WriteOffsetFile()
{ 
  for( unsigned int id = m_ZTileStart; id <= m_ZTileEnd; id++  )
  {    
    std::stringstream filename;
    filename << m_OffsetFilePath << id+1 << ".txt";
    std::ofstream os ( filename.str().c_str() );

    if ( !os )
    {
      std::cout << "error in offset file opening" << std::endl;
      return;
    }

    for( unsigned int i = 0; i < ImageDimension; i++ )
    {
      for( unsigned int j = 0; j < m_TileNumber[2]; j++ )
      {
        if ( j == id+1 )
        {
          os << m_SharedData->m_TileOffset[i][j] << ' ';
        }
        else
        {
          os << 0.0 << ' ';
        }
      }
      os << std::endl;
    }

    os.close();
  }
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

  std::vector< std::string > relevantFilenameVector;
  searchStringCH << m_ChannelPrefix << m_ChannelNumber;
  searchStringXYZT << std::setfill( '0' ) << std::setw( 4 ) << m_TimePoint << "t";
  for ( unsigned int m = 0; m < directory->GetNumberOfFiles(); m++)
  {
    filename = directory->GetFile( m );

    if ( ( filename.find( searchStringCH.str() ) != std::string::npos ) &&
         ( filename.find( searchStringXYZT.str() ) != std::string::npos ) )
    {
      relevantFilenameVector.push_back( filename );
    }
  }

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
        searchStringXYZT.str( std::string() );
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << i << "x_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << j << "y_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 3 ) << k << "z_";
        searchStringXYZT << std::setfill( '0' ) << std::setw( 4 ) << m_TimePoint << "t";

        for ( unsigned int m = 0; m < relevantFilenameVector.size(); m++)
        {
          filename = relevantFilenameVector[m];

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

//    const DictionaryType & dictionary = reader->GetImageIO()->GetMetaDataDictionary();
//    typename DictionaryType::ConstIterator itr = dictionary.Begin();
//    typename DictionaryType::ConstIterator end = dictionary.End();

//    unsigned int count = 0;
//    while( itr != end )
//    {
//      std::cout << count++ << std::endl;
//      typename MetaDataObjectBase::Pointer entry = itr->second;
//      typename MetaDataStringType::Pointer entryValue = dynamic_cast<MetaDataStringType *> (
//            entry.GetPointer() );

//      if ( entryValue )
//      {
//        std::string tagkey = itr->first;
//        std::cout << "MetaDataDict " << tagkey << ": " << entryValue->GetMetaDataObjectValue() << std::endl;
//      }
//      ++itr;
//    }

    m_TileDimension = currentImage->GetLargestPossibleRegion().GetSize();

    m_TileSpacing[0] = m_TileSize[0]/m_TileDimension[1];
    m_TileSpacing[1] = m_TileSize[1]/m_TileDimension[0];
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
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_SharedData->m_TileOffset[i].resize( m_TileNumber[2], 0.0 );
    m_SharedData->m_TileEffectiveOffset[i].resize( m_TileNumber[2], 0.0 );
  }

  // Read offset file
  if ( !m_OffsetFilePath.empty() )
  {
    ReadOffsetFile();
  }

  if ( m_RegisterZTiles )
  {
    for( unsigned int i = 0; i < m_SearchRadius.size(); i++ )
    {
      RegisterTiles( m_SearchRadius[i], m_StepLength[i] );
    }

    // Write out the offsets
    std::cout << "Offset file path: " << m_OffsetFilePath << std::endl;
    if ( !m_OffsetFilePath.empty() )
    {
      WriteOffsetFile();
    }
  }

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
