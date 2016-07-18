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

#ifndef __itkFillROIImageFilter_txx
#define __itkFillROIImageFilter_txx

#include "itkObjectFactory.h"
#include "itkFillROIImageFilter.h"

namespace itk
{

template< class TInputImage >
FillROIImageFilter< TInputImage >
::FillROIImageFilter()
{
  m_ZTile = 0;
  m_SingleZFill = false;
  this->SetNumberOfRequiredInputs(1);
  Superclass::GenerateInputRequestedRegion();

  ImagePointer input = const_cast< ImageType * >( this->GetInput() );
  if ( input )
  {
    input->SetRequestedRegion( input->GetLargestPossibleRegion() );
  }
}

template< class TInputImage >
void
FillROIImageFilter< TInputImage >::
OverlapRegion( ImagePointer A, ImagePointer B,
  RegionType& rA, RegionType& rB )
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


template< class TInputImage >
void
FillROIImageFilter< TInputImage >::
BeforeThreadedGenerateData()
{
  this->GenerateOutputInformation();
  this->AllocateOutputs();
  this->ReleaseDataFlagOn();

  const ImageType *outputPtr = this->GetOutput();
  const ImageRegionSplitterBase * splitter = this->GetImageRegionSplitter();
  m_NumOfValidThreads = splitter->GetNumberOfSplits( outputPtr->GetRequestedRegion(), this->GetNumberOfThreads() );

  ImagePointer m_ROIImage = this->GetOutput();
  PointType m_ROIOrigin = m_ROIImage->GetOrigin();
  SpacingType m_TileSpacing = m_ROIImage->GetSpacing();
  RegionType m_ROI = m_ROIImage->GetLargestPossibleRegion();
  
  unsigned int m_TileNumber[3];
  for( unsigned int k = 0; k < ImageDimension; k++ )
  {
    m_ScanStart[k] = 0;
    m_TileNumber[k] = m_SharedData->m_TileCover[k][0][0].size();
    m_ScanEnd[k] = m_TileNumber[k] - 1;
  }
  std::cout << "Allocated ROI image" << std::endl;

  if ( m_SingleZFill )
  {
    m_ScanStart[2] = m_ScanEnd[2] = m_ZTile;
    return;
  }
  
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
      //std::cout << k << ' ' << i << ' ' << m_SharedData->m_TileCoverStart[k][i] << ' '
      //          << m_SharedData->m_TileCoverEnd[k][i] << std::endl;
      if ( ( beginCorner >= m_SharedData->m_TileCover[k][0][0][i] - 0.0001 ) &&
           ( beginCorner <= m_SharedData->m_TileCover[k][1][0][i] + 0.0001 ) &&
           ( scanStartVal >= m_SharedData->m_TileCover[k][0][0][i] ) )
      {
        m_ScanStart[k] = i;
        scanStartVal =  m_SharedData->m_TileCover[k][0][0][i];
      }

      if ( ( endCorner >= m_SharedData->m_TileCover[k][0][0][i] - 0.001 ) &&
           ( endCorner <= m_SharedData->m_TileCover[k][1][0][i] + 0.001 ) &&
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

    std::cout << m_ScanStart[k] << ' ' << m_ScanEnd[k] << std::endl;
  }
}


template< class TInputImage >
void
FillROIImageFilter< TInputImage >::
ThreadedGenerateData(const RegionType &windowRegion, ThreadIdType threadId)
{
  ImagePointer m_ROIImage = this->GetOutput();
  SpacingType m_TileSpacing = m_ROIImage->GetSpacing();
  RegionType m_ROI = m_ROIImage->GetLargestPossibleRegion();

  // Start a loop that will read all the tiles from zScanStart to zScanEnd
  PointType currentTileOrigin;
  RegionType currentTileRegion, currentTileOverlapRegion, roiOverlapRegion;
  RegionType roi;
  SizeType clipTileSize;
  IndexType clipTileIndex;
  PointType clipTileOrigin;

  unsigned int counter = 0;

  for( unsigned int i = m_ScanStart[0]; i <= m_ScanEnd[0]; i++ )
  {
    for( unsigned int j = m_ScanStart[1]; j <= m_ScanEnd[1]; j++ )
    {
      for( unsigned int k = m_ScanStart[2]; k <= m_ScanEnd[2]; k++, counter++ )
      {
        //std::cout << i << ' ' << j << ' ' << k << std::endl;
        if ( counter%(m_NumOfValidThreads) == threadId )
        {
          //std::cout << counter << ' ' << counter%(m_NumOfValidThreads) << std::endl;
          std::string filename = m_SharedData->m_TileFileNameArray[i][j][k];
          //std::cout << filename.c_str() << std::endl;
          if  ( ! filename.empty() )
          {
            currentTileOrigin[0] = m_SharedData->m_TileCover[0][0][0][i];
            currentTileOrigin[1] = m_SharedData->m_TileCover[1][0][0][j];
            currentTileOrigin[2] = m_SharedData->m_TileCover[2][0][0][k];

            //std::cout << "Current Tile Origin " << currentTileOrigin << std::endl;

            clipTileOrigin[0] = m_SharedData->m_TileCover[0][0][1][i];
            clipTileOrigin[1] = m_SharedData->m_TileCover[1][0][1][j];
            clipTileOrigin[2] = m_SharedData->m_TileCover[2][0][1][k];

            //std::cout << "Clip Tile Origin " << clipTileOrigin << std::endl;

            clipTileSize[0] = 1 + static_cast<SizeValueType>(
                          ( m_SharedData->m_TileCover[0][1][1][i] - m_SharedData->m_TileCover[0][0][1][i] )/m_TileSpacing[0] );
            clipTileSize[1] = 1 + static_cast<SizeValueType>(
                          ( m_SharedData->m_TileCover[1][1][1][j] - m_SharedData->m_TileCover[1][0][1][j] )/m_TileSpacing[1] );
            clipTileSize[2] = 1 + static_cast<SizeValueType>(
                          ( m_SharedData->m_TileCover[2][1][1][k] - m_SharedData->m_TileCover[2][0][1][k] )/m_TileSpacing[2] );

            //std::cout << "Clip Tile Size " << clipTileSize << std::endl;

            ImagePointer tileImage = ExtractCorrectedAndFlippedTile( filename );
            tileImage->SetOrigin( currentTileOrigin );
            SizeType m_TileDimension = tileImage->GetLargestPossibleRegion().GetSize();

            //std::cout << "Extraction complete" << std::endl;

            tileImage->TransformPhysicalPointToIndex( clipTileOrigin, clipTileIndex );

            for( unsigned int m = 0; m < ImageDimension; m++ )
            {
              if ( clipTileIndex[m] + clipTileSize[m] > m_TileDimension[m] )
              {
                clipTileSize[m] = m_TileDimension[m] - clipTileIndex[m] - 1;
              }
            }

            roi.SetSize( clipTileSize );
            roi.SetIndex( clipTileIndex );


            //std::cout << "Tile region:" << tileImage->GetLargestPossibleRegion() << std::endl;
            //std::cout << "Clip tile roi: " << roi << std::endl;

            // Extract ROI
            ImagePointer currentImage = tileImage;
            if ( !m_SingleZFill )
            {
              ROIFilter3DPointer roiFilter = ROIFilter3DType::New();
              roiFilter->SetRegionOfInterest( roi );
              roiFilter->SetInput( tileImage );
              roiFilter->Update();
              currentImage = roiFilter->GetOutput();
              currentImage->DisconnectPipeline();
            }
            else
            {
	      //std::cout << m_SingleZFill << ' ' << m_ZTile << std::endl;
	    }
            currentTileRegion = currentImage->GetLargestPossibleRegion();

            //std::cout << "ROI filtering complete " << std::endl;

            OverlapRegion( currentImage , m_ROIImage,
                           currentTileOverlapRegion, roiOverlapRegion );

            //std::cout << "ROI Image origin: " << m_ROIImage->GetOrigin() << std::endl;
            //std::cout << "ROI extent: " << m_ROI << std::endl;
            //std::cout << "ROI Image region: " << roiOverlapRegion << std::endl;

            //std::cout << "Tile origin: " << currentImage->GetOrigin() << std::endl;
            //std::cout << "Current Tile extent: " << currentImage->GetLargestPossibleRegion() << std::endl;
            //std::cout << "Current Tile Region: " << currentTileOverlapRegion << std::endl;

            // Using these images, fill up roiImage

            if ( m_ROI.IsInside( roiOverlapRegion ) &&
                 currentTileRegion.IsInside( currentTileOverlapRegion ) )
            {// Clipping may eliminate overlap
              IteratorType rIt( m_ROIImage, roiOverlapRegion );
              rIt.GoToBegin();
              IteratorType tIt( currentImage, currentTileOverlapRegion );
              tIt.GoToBegin();

              PixelType p;
              while( !tIt.IsAtEnd() )
              {
                rIt.Set( tIt.Get() );//rIt.Get() +
                ++tIt;
                ++rIt;
              }
            }
          }
        }
      }
    }
  }
}


template< class TInputImage >
typename FillROIImageFilter< TInputImage >::ImagePointer
FillROIImageFilter< TInputImage >::
ExtractCorrectedAndFlippedTile( std::string& filename )
{
  ImagePointer m_ROIImage = this->GetOutput();
  PointType m_ROIOrigin = m_ROIImage->GetOrigin();
  SpacingType m_TileSpacing = m_ROIImage->GetSpacing();
  RegionType m_ROI = m_ROIImage->GetLargestPossibleRegion();

  ReaderPointer m_Reader = ReaderType::New();
  m_Reader->SetFileName( filename.c_str() );
  m_Reader->Update();

  ImagePointer cImage = m_Reader->GetOutput();
  cImage->DisconnectPipeline();

  PixelType p;
  double num, den;
  if ( m_SharedData->m_CorrectionImage )
  {
    //std::cout << "Correction used" << std::endl;
    IteratorType cIt( cImage, cImage->GetLargestPossibleRegion() );
    cIt.GoToBegin();
    RIteratorType corrIt( m_SharedData->m_CorrectionImage,
                          m_SharedData->m_CorrectionImage->GetLargestPossibleRegion() );
    corrIt.GoToBegin();
    while( !cIt.IsAtEnd() )
    {
      if ( corrIt.IsAtEnd() )
      {
        corrIt.GoToBegin();
      }

      p = static_cast<PixelType>( m_SharedData->m_CorrectionThreshold );
      if ( cIt.Get() > m_SharedData->m_CorrectionThreshold )
      {
        num = double(cIt.Get()) - m_SharedData->m_CorrectionThreshold;
        den = corrIt.Get();
        if ( den < 0.01 )
        {
          den = 0.01;
        }
        p += static_cast<PixelType>( num/den );
        cIt.Set( p );
      }
      ++cIt;
      ++corrIt;
    }
  }

  FixedArray<unsigned int, 3> axesOrder;
  axesOrder[0] = 1;
  axesOrder[1] = 0;
  axesOrder[2] = 2;

  PermuteAxesFilterPointer pAFilter = PermuteAxesFilterType::New();
  pAFilter->SetInput( cImage );
  pAFilter->SetOrder( axesOrder );
  pAFilter->Update();
  ImagePointer pImage = pAFilter->GetOutput();

  ImagePointer currenTInputImage = ImageType::New();
  currenTInputImage->SetOrigin( pImage->GetOrigin() );
  currenTInputImage->SetSpacing( m_TileSpacing );
  currenTInputImage->SetRegions( pImage->GetLargestPossibleRegion() );
  currenTInputImage->Allocate();
  currenTInputImage->FillBuffer( 0 );

  PixelType q;
  IteratorType pIt( pImage, pImage->GetLargestPossibleRegion() );
  IteratorType currentIt( currenTInputImage,
                          currenTInputImage->GetLargestPossibleRegion() );
  while( !pIt.IsAtEnd() )
  {
    currentIt.Set( pIt.Get() );
    ++pIt;
    ++currentIt;
  }

  ImagePointer deconvImage;
  if ( m_SharedData->m_PSF )
  {
    //std::cout << "Deconvolution used" << std::endl;
    DeconvolutionFilterPointer deconvolutionFilter = DeconvolutionFilterType::New();
    deconvolutionFilter->SetInput( currenTInputImage );
    deconvolutionFilter->SetKernelImage( m_SharedData->m_PSF );
    deconvolutionFilter->NormalizeOn();
    deconvolutionFilter->SetNumberOfIterations( m_SharedData->m_DeconvolutionIterations );
    deconvolutionFilter->Update();
    deconvImage = deconvolutionFilter->GetOutput();
    deconvImage->DisconnectPipeline();
  }
  else
  {
    deconvImage = currenTInputImage;
  }

  return deconvImage;
}


/** Print Self information */
template<class TInputImage >
void
FillROIImageFilter< TInputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace

#endif
