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

#ifndef __itkExtractStitchTransformImageFilter_txx
#define __itkExtractStitchTransformImageFilter_txx

#include "itkExtractStitchTransformImageFilter.h"

static itk::SimpleFastMutexLock m_MutexESTIF;

namespace itk
{
template < class TInputImage >
ExtractStitchTransformImageFilter< TInputImage >
::ExtractStitchTransformImageFilter()
{
  m_OffsetFilePath = "";
  m_SharedData = ITK_NULLPTR;
  m_ZTileStart = 0;
  m_ZTileEnd = 0;
  m_NumberOfThreads = 1;
}


template < class TInputImage >
void
ExtractStitchTransformImageFilter< TInputImage >::
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


template < class TInputImage >
void
ExtractStitchTransformImageFilter< TInputImage >::
BeginRegister()
{
  // Read offset file
  if ( !m_OffsetFilePath.empty() )
  {
    ReadOffsetFile();
  }

  if ( m_ZTileEnd > m_SharedData->m_TileNumber[2]-2 )
  {
    m_ZTileEnd = m_SharedData->m_TileNumber[2]-2;
  }

  ThreadStruct str;
  str.Filter  = this;

  ThreaderPointer threader = ThreaderType::New();
  threader->SetNumberOfThreads( m_NumberOfThreads );
  threader->SetSingleMethod( this->ThreaderCallback, &str );
  threader->SingleMethodExecute();

  // Write out the offsets
  if ( !m_OffsetFilePath.empty() )
  {
    WriteOffsetFile();
  }
}


template < class TInputImage >
ITK_THREAD_RETURN_TYPE
ExtractStitchTransformImageFilter< TInputImage >::
ThreaderCallback(void * arg)
{
  unsigned int ThreadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;

  ThreadStruct * str =
  (ThreadStruct *) (((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  m_MutexESTIF.Lock();
  std::cout << ThreadId << " == " << "Start processing" << std::endl;
  m_MutexESTIF.Unlock();

  for( unsigned int i = 0; i < str->Filter->m_SearchRadius.size(); i++ )
  {
    for( unsigned int ztile = str->Filter->m_ZTileStart; ztile <= str->Filter->m_ZTileEnd; ztile++ )
    {
      unsigned int rem = ( ztile ) % ( str->Filter->m_NumberOfThreads );
      if ( rem == ThreadId )
      {
        str->Filter->RegisterTiles( str->Filter->m_SearchRadius[i], str->Filter->m_StepLength[i], ztile );
      }
    }
  }

  m_MutexESTIF.Lock();
  std::cout << ThreadId << " == " << "End processing" << std::endl;
  m_MutexESTIF.Unlock();

  return ITK_THREAD_RETURN_VALUE;
}


template < class TInputImage >
void
ExtractStitchTransformImageFilter< TInputImage >::
RegisterTiles( float searchRadius, float stepLength, unsigned int ztile )
{
  unsigned int i_,j_;
  {
    bool tileForRegistration = false;
    i_ = m_SharedData->m_TileNumber[0]/2;
    j_ = m_SharedData->m_TileNumber[1]/2;
    if ( ( !m_SharedData->m_TileFileNameArray[i_][j_][ztile].empty() ) &&
    ( !m_SharedData->m_TileFileNameArray[i_][j_][ztile+1].empty() ) )
    {
      tileForRegistration = true;
    }

    for( unsigned int i = 0; i < m_SharedData->m_TileNumber[0]; i++ )
    {
      for( unsigned int j = 0; j < m_SharedData->m_TileNumber[1]; j++ )
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
    morigin[2] = m_SharedData->m_TileSize[2] - m_SharedData->m_TileOverlap[2]
        + m_SharedData->m_TileOffset[2][ztile+1];

    staticImage->SetOrigin( sorigin );
    movingImage->SetOrigin( morigin );

    staticImage->SetSpacing( m_SharedData->m_TileSpacing );
    movingImage->SetSpacing( m_SharedData->m_TileSpacing );
    
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
            //std::cout << "*****" << besti<< ' ' << bestj << ' '
            //          << bestk << ' ' << bestValue << std::endl;
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


template < class TInputImage >
void
ExtractStitchTransformImageFilter< TInputImage >::
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


template < class TInputImage >
void
ExtractStitchTransformImageFilter< TInputImage >::
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
      for( unsigned int j = 0; j < m_SharedData->m_TileNumber[2]; j++ )
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

} /* end namespace itk */

#endif
