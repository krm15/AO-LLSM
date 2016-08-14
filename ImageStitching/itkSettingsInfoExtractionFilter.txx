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
  m_SettingFieldName.resize( 100, "" );
  m_SettingFieldValue.resize( 100, "" );

  m_NumberOfTiles = 1;
  m_StitchedImage = ITK_NULLPTR;
  m_ChannelPrefix = "_ch";
  m_OffsetFilePath = "";
  m_SharedData = ITK_NULLPTR;
  m_SampleScan = false;

  tileAxesOrder[0] = 0;
  tileAxesOrder[1] = 1;
  tileAxesOrder[2] = 2;
}


template < class TValueType, class TInputImage >
void
SettingsInfoExtractionFilter< TValueType, TInputImage >::
ReadOffsetFile()
{
  for( unsigned int id = 1; id < m_SharedData->m_TileNumber[2]; id++  )
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
      for( unsigned int j = 0; j < m_SharedData->m_TileNumber[2]; j++ )
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
    for( unsigned int j = 1; j < m_SharedData->m_TileNumber[2]; j++ )
    {
      m_SharedData->m_TileEffectiveOffset[i][j] =
          m_SharedData->m_TileEffectiveOffset[i][j-1] +
          m_SharedData->m_TileOffset[i][j];
    }
  }

  // Flip effective offsets
  double temp;
  for( unsigned int j = 1; j < m_SharedData->m_TileNumber[2]; j++ )
  {
    temp = m_SharedData->m_TileEffectiveOffset[0][j];
    m_SharedData->m_TileEffectiveOffset[0][j] = m_SharedData->m_TileEffectiveOffset[1][j];
    m_SharedData->m_TileEffectiveOffset[1][j] = temp;
  }
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
        m_SharedData->m_TileCover[i][j][k].resize( m_SharedData->m_TileNumber[i] );
      }
    }
  }

  for( unsigned int j = 0; j < ImageDimension; j++ )
  {
    for( unsigned int i = 0; i < m_NumberOfTiles; i++ )
    {
      unsigned int jj = tileAxesOrder[j];
      unsigned int temp =  m_TileInfoValue[i][jj];// range is m_TileNumber[j]
      m_SharedData->m_TileCover[j][0][0][temp] = m_TransformedTileInfoValue[i][jj+3];
      m_SharedData->m_TileCover[j][1][0][temp] = m_TransformedTileInfoValue[i][jj+6];
    }
  }

  // HACK
  for( unsigned int j = 0; j < ImageDimension; j++ )
  {
    bool sign = true;
    if ( !m_SampleScan )
    {
      if ( m_SharedData->m_TileCover[j][0][0][0] > m_SharedData->m_TileCover[j][0][0][1] )
      {
        sign = false;
      }
    }

    for( unsigned int k = 0; k < m_SharedData->m_TileNumber[j]; k++ )
    {
       m_SharedData->m_TileCover[j][0][0][k] = ( k * ( m_SharedData->m_TileSize[j] -  m_SharedData->m_TileOverlap[j]) );

       if (!sign)
       {
         m_SharedData->m_TileCover[j][0][0][k] = -m_SharedData->m_TileCover[j][0][0][k];
       }

       m_SharedData->m_TileCover[j][1][0][k] = m_SharedData->m_TileCover[j][0][0][k] + m_SharedData->m_TileSize[j];
    }
  }


  for( unsigned int j = 0; j < ImageDimension; j++ )
  {
    m_MinimumStart[j] = 1000000.0;
    m_MaximumEnd[j]   = -m_MinimumStart[j];
    for( unsigned int k = 0; k < m_SharedData->m_TileNumber[j]; k++ )
    {
      // Clipping should be based on offsets as well
      m_SharedData->m_TileCover[j][0][1][k] = m_SharedData->m_TileCover[j][0][0][k]
          + 0.5*m_SharedData->m_TileOverlap[j];
      m_SharedData->m_TileCover[j][1][1][k] = m_SharedData->m_TileCover[j][1][0][k]
          - 0.5*m_SharedData->m_TileOverlap[j];

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
    unsigned int firstTile(0), lastTile(m_SharedData->m_TileNumber[j]-1);

    for( unsigned int k = 0; k < m_SharedData->m_TileNumber[j]; k++ )
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
ReadTileInfo( std::istream& os )
{
  for( unsigned int i = 0; i < m_SettingFieldName.size(); i++ )
  {
    if ( m_SettingFieldName[i] == "# Subvolume X")
    {
        m_SharedData->m_TileNumber[0] = atoi( m_SettingFieldValue[i].c_str() );
        //std::cout << "Tile Number[0] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "# Subvolume Y")
    {
        m_SharedData->m_TileNumber[1] = atoi( m_SettingFieldValue[i].c_str() );
        //std::cout << "Tile Number[1] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "# Subvolume Z")
    {
        m_SharedData->m_TileNumber[2] = atoi( m_SettingFieldValue[i].c_str() );
        //std::cout << "Tile Number[2] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "Size of Z Stack X (um)")
    {
        m_SharedData->m_TileSize[0] = atoi( m_SettingFieldValue[i].c_str() );
        //std::cout << "Tile size[0] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "Size of Z Stack Y (um)")
    {
        m_SharedData->m_TileSize[1] = atoi( m_SettingFieldValue[i].c_str() );
        //std::cout << "Tile size[1] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "Size of Z Stack Z (um)")
    {
        m_SharedData->m_TileSize[2] = atoi( m_SettingFieldValue[i].c_str() );
        //std::cout << "Tile size[2] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "Overlap X (um)")
    {
        m_SharedData->m_TileOverlap[0] = atof( m_SettingFieldValue[i].c_str() );
        //std::cout << "Tile overlap[0] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "Overlap Y (um)")
    {
        m_SharedData->m_TileOverlap[1] = atof( m_SettingFieldValue[i].c_str() );
        //std::cout << "Tile overlap[1] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "Overlap Z (um)")
    {
        m_SharedData->m_TileOverlap[2] = atof( m_SettingFieldValue[i].c_str() );
        //std::cout << "Tile overlap[2] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "Z Motion")
    {
      if ( m_SettingFieldValue[i] == "Sample piezo" )
      {
        m_SampleScan = true;
      }
      //std::cout << "Tile overlap[2] " << m_SettingFieldValue[i] << std::endl;
    }
    if ( m_SettingFieldName[i] == "PC name")
    {
      m_ScopeName = m_SettingFieldValue[i];
    }
  }

  // Flip m_TileNumber and m_TileOverlap for sample scan
  if ( m_SampleScan )
  {
    tileAxesOrder[0] = 1;
    tileAxesOrder[1] = 0;

    unsigned int temp1 = m_SharedData->m_TileNumber[0];
    m_SharedData->m_TileNumber[0] = m_SharedData->m_TileNumber[1];
    m_SharedData->m_TileNumber[1] = temp1;

    double temp2 = m_SharedData->m_TileOverlap[0];
    m_SharedData->m_TileOverlap[0] = m_SharedData->m_TileOverlap[1];
    m_SharedData->m_TileOverlap[1] = temp2;
  }

  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_NumberOfTiles *= m_SharedData->m_TileNumber[i];
  }

  std::string value, line;
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
      m_TransformedTileInfoValue[i][j+ImageDimension] = newCenter[j] - 0.5*m_SharedData->m_TileSize[j];
      m_TransformedTileInfoValue[i][j+2*ImageDimension] = newCenter[j] + 0.5*m_SharedData->m_TileSize[j];
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

  searchStringCH << m_ChannelPrefix << m_ChannelNumber;
  searchStringXYZT << std::setfill( '0' ) << std::setw( 4 ) << m_TimePoint << "t";

  m_SharedData->m_TileFileNameArray.resize( m_SharedData->m_TileNumber[0] );
  for( unsigned int i = 0; i < m_SharedData->m_TileNumber[0]; i++ )
  {
    m_SharedData->m_TileFileNameArray[i].resize( m_SharedData->m_TileNumber[1] );
    for( unsigned int j = 0; j < m_SharedData->m_TileNumber[1]; j++ )
    {
      m_SharedData->m_TileFileNameArray[i][j].resize( m_SharedData->m_TileNumber[2] );
      for( unsigned int k = 0; k < m_SharedData->m_TileNumber[2]; k++ )
      {
        m_SharedData->m_TileFileNameArray[i][j][k] = std::string();
      }
    }
  }

  unsigned int xp, yp, zp, pos, p;
  unsigned int totalPixelCount = 0;
  unsigned int m_TrueCountOfTiles = 0;
  std::vector< unsigned int > histogram;
  IndexType index, index2;
  histogram.resize(300000, 0);
  for ( unsigned int m = 0; m < directory->GetNumberOfFiles(); m++)
  {
    //std::cout << "m: " << m << std::endl;
    filename = directory->GetFile( m );

    if ( ( filename.find( searchStringCH.str() ) != std::string::npos ) &&
         ( filename.find( searchStringXYZT.str() ) != std::string::npos ) )
    {
      pos = filename.find_last_of("x");
      index2[0] = atoi( filename.substr( pos-3, 3 ).c_str() );
      pos = filename.find_last_of("y");
      index2[1] = atoi( filename.substr( pos-3, 3 ).c_str() );
      pos = filename.find_last_of("z");
      index2[2] = atoi( filename.substr( pos-3, 3 ).c_str() );

      if ( m_SampleScan )
      {
        index[0] = index2[1];
        index[1] = m_SharedData->m_TileNumber[0] - index2[0];
        index[2] = m_SharedData->m_TileNumber[2] - index2[2] - 1;
      }
      else
      {
        index = index2;
      }

      //std::cout << index << ' ' << index2 << std::endl;

      m_SharedData->m_TileFileNameArray[index[0]][index[1]][index[2]] = m_TileDirectory + filename;

      unsigned int lastindex = filename.find_last_of("d");
      std::string filename_mip = m_TileDirectory + "MIPs/" + filename.substr(0, lastindex) + "MIP_z.tif";

      std::string filename_raw = m_SettingsDirectory + filename.substr(0, lastindex-1) + ".tif";

      //std::ifstream infile( filename_mip.c_str() );
      //if ( infile )
      {
        //infile.close();
        RReaderPointer reader = RReaderType::New();
        reader->SetFileName ( filename_mip.c_str() );
        reader->Update();

        RImagePointer img = reader->GetOutput();
        totalPixelCount += img->GetLargestPossibleRegion().GetNumberOfPixels();

        RIteratorType It( img, img->GetLargestPossibleRegion() );
        It.GoToBegin();
        while( !It.IsAtEnd() )
        {
          p = static_cast<unsigned int>(It.Get());
          histogram[p]++;
          ++It;
        }
      }

      // Read the associated tags of this filename
      TIFF* image = TIFFOpen(filename_raw.c_str(), "r");
      float *cenx, *ceny, *cenz;
      unsigned int count;
      TIFFGetField(image, 40000, &count, &cenx);
      TIFFGetField(image, 40001, &count, &ceny);
      TIFFGetField(image, 40002, &count, &cenz);
      //std::cout << filename_raw.c_str() << ' ' << *cenx << ' ' << *ceny << ' ' << *cenz << std::endl;

      if ( !ChannelNameSet )
      {
        pos = filename.find( searchString );
        m_ChannelName = filename.substr( pos-3, 5 );
        ChannelNameSet = true;
        m_SampleName = m_SharedData->m_TileFileNameArray[index[0]][index[1]][index[2]];
      }

      m_TrueCountOfTiles++;
    }
  }
  m_NumberOfTiles = m_TrueCountOfTiles;

  unsigned int topPercentOfPixels = 0.03 * totalPixelCount;
  unsigned int max = histogram.size() - 1;
  unsigned int cumsum = 0;
  while( ( max > 0 ) && ( cumsum < topPercentOfPixels ) )
  {
      cumsum += histogram[max];
      max--;
  }

  // Do some analysis here to find a good m

  if ( max > 0 )
  {
    m_SharedData->m_ScalingFactor = 65535/max;
  }
  else
  {
    m_SharedData->m_ScalingFactor = 1.0;
  }
}


template < class TValueType, class TInputImage >
std::istream&
SettingsInfoExtractionFilter< TValueType, TInputImage >::
safeGetline(std::istream& is, std::string& t)
{
    t.clear();

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(;;) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            if(t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t += (char)c;
        }
    }
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
  std::string value1, value2, line1, line2;
  safeGetline ( os, line1 );
  std::stringstream nameStream( line1 );

  //std::cout << "Line1: " << line1 << std::endl;

  safeGetline ( os, line2);
  std::stringstream valueStream( line2 );

  //std::cout << "Line2: " << line2 << std::endl;

  // First two lines is 29-34 fields of data
  unsigned int i = 0;
  while( std::getline ( nameStream, value1, ',' ) )
  {
    m_SettingFieldName[i] = value1;

    //std::cout << value1 << ' ';
    std::getline ( valueStream, value2, ',' );
    m_SettingFieldValue[i] = value2;
    //std::cout << value2 << std::endl;
    i++;
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

  m_SharedData->m_TileDimension = currentImage->GetLargestPossibleRegion().GetSize();

  // Store m_TileDimension correctly since image is flipped
  unsigned int temp = m_SharedData->m_TileDimension[0];
  m_SharedData->m_TileDimension[0] = m_SharedData->m_TileDimension[1];
  m_SharedData->m_TileDimension[1] = temp;

  m_SharedData->m_TileSpacing[0] = m_SharedData->m_TileSize[0]/m_SharedData->m_TileDimension[0];
  m_SharedData->m_TileSpacing[1] = m_SharedData->m_TileSize[1]/m_SharedData->m_TileDimension[1];
  m_SharedData->m_TileSpacing[2] = m_SharedData->m_TileSize[2]/m_SharedData->m_TileDimension[2];
  std::cout << "Read tile dimensions" << std::endl;

  // Create a stitched image
  CreateStitchedImage();
  std::cout << "Estimate stitched image specification" << std::endl;

  // Register files
  for( unsigned int i = 0; i < ImageDimension; i++ )
  {
    m_SharedData->m_TileOffset[i].resize( m_SharedData->m_TileNumber[2], 0.0 );
    m_SharedData->m_TileEffectiveOffset[i].resize( m_SharedData->m_TileNumber[2], 0.0 );
  }

  // Read offset file
  if ( !m_OffsetFilePath.empty() )
  {
    ReadOffsetFile();
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
    m_StitchDimension[i] = m_StitchSize[i]/m_SharedData->m_TileSpacing[i];
  }

  m_StitchedImage = ImageType::New();
  m_StitchedImage->SetOrigin( m_StitchOrigin );
  m_StitchedImage->SetSpacing( m_SharedData->m_TileSpacing );
  m_StitchedImage->SetRegions( m_StitchRegion );
}

} /* end namespace itk */

#endif
