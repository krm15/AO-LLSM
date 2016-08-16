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

#ifndef __itkStitchingSharedData_txx
#define __itkStitchingSharedData_txx

#include "itkObjectFactory.h"
#include "itkStitchingSharedData.h"

namespace itk
{

template< class TInputImage >
StitchingSharedData< TInputImage >
::StitchingSharedData()
{
  m_CorrectionImage = ITK_NULLPTR;
  m_CorrectionThreshold = 1200;
  m_CorrectionVariance = 2.0;
  m_CorrectionFilename = "";
  m_PSF = ITK_NULLPTR;
  m_PSFPath = "";
  m_DeconvolutionIterations = 15;
  m_ScalingFactor = 1.0;
}


template < class TInputImage >
void
StitchingSharedData< TInputImage >::
ComputeScalingFactor( std::string& iDirName )
{
  DirectoryPointer directory = DirectoryType::New();
  directory->Load( iDirName.c_str() );

  std::string filename, filename2;
  std::string searchString = "_MIP_z.tif";

  unsigned int histogramSize = 150000000;
  std::vector< unsigned int > histogram;
  histogram.resize( histogramSize, 0 );
  unsigned int m_TrueCountOfTiles = 0;

  unsigned int totalPixelCount = 0;
  unsigned int numOfFiles = directory->GetNumberOfFiles();

  unsigned int p;
  double _max = 0.0;
  for ( unsigned int m = 0; m < numOfFiles; m++)
  {
    //std::cout << "m: " << m << std::endl;
    filename = directory->GetFile( m );
    if ( filename.find( searchString ) != std::string::npos )
    {
      filename2 = iDirName + filename;
      std::ifstream infile( filename2.c_str() );
      if ( infile )
      {
        infile.close();
        //std::cout << filename2 << std::endl;

        RReaderPointer reader = RReaderType::New();
        reader->SetFileName ( filename2.c_str() );
        reader->Update();

        RImagePointer img = reader->GetOutput();
        totalPixelCount += img->GetLargestPossibleRegion().GetNumberOfPixels();

        RIteratorType It( img, img->GetLargestPossibleRegion() );
        It.GoToBegin();
        while( !It.IsAtEnd() )
        {
//          if ( It.Get() > _max )
//          {
//            _max =  It.Get();
//          }

          p = static_cast<unsigned int>( It.Get() );
          if ( p >= histogramSize )
          {
            p = histogramSize-1;
          }
          histogram[p]++;
          ++It;
        }
        m_TrueCountOfTiles++;
      }
    }
  }

  //std::cout << "Max pixel intensity: " << _max << std::endl;

  if ( totalPixelCount > 0 )
  {
    unsigned int topPercentOfPixels = 0.03 * totalPixelCount;
    unsigned int maxIndex = histogramSize - 1;
    unsigned int cumsum = 0;
    while( ( maxIndex > 0 ) && ( cumsum < topPercentOfPixels ) )
    {
      cumsum += histogram[maxIndex];
      maxIndex--;
    }

    if ( maxIndex > 0 )
    {
     m_ScalingFactor = double(65535)/( (double)maxIndex );
    }
  }
  std::cout << "Scaling Factor: " << m_ScalingFactor << std::endl;
}


template < class TInputImage >
void
StitchingSharedData< TInputImage >::
ReadPSFImage()
{
  if ( m_PSFPath.empty() )
  {
    return;
  }

  // Read the correction image
  ReaderPointer reader = ReaderType::New();
  reader->SetFileName ( m_PSFPath );
  reader->Update();

  m_PSF = reader->GetOutput();
  m_PSF->DisconnectPipeline();

  SpacingType sp;
  sp[0] = m_TileSpacing[0];
  sp[1] = m_TileSpacing[1];
  sp[2] = m_TileSpacing[2];
  m_PSF->SetSpacing( sp );
}


template < class TInputImage >
void
StitchingSharedData< TInputImage >::
ReadCorrectionImage()
{
  RSpacingType sp;
  sp[0] = m_TileSpacing[0];
  sp[1] = m_TileSpacing[1];

  // Read the correction image
  RReaderPointer reader = RReaderType::New();
  reader->SetFileName ( m_CorrectionFilename.c_str() );
  reader->Update();

  RImagePointer inputImage = reader->GetOutput();
  inputImage->DisconnectPipeline();
  inputImage->SetSpacing( sp );

  GaussianFilterPointer gaussianFilter = GaussianFilterType::New();
  gaussianFilter->SetInput( inputImage );
  gaussianFilter->SetVariance( m_CorrectionVariance );
  gaussianFilter->SetUseImageSpacingOn();
  gaussianFilter->Update();

  RImagePointer currentImage = gaussianFilter->GetOutput();
  currentImage->DisconnectPipeline();
  currentImage->SetSpacing( sp );

  if ( !currentImage )
  {
      return;
  }

  RRegionType rroi;

  RSizeType rsize = currentImage->GetLargestPossibleRegion().GetSize();

  std::cout << "Correction map size: " << rsize << std::endl;
  //std::cout << m_TileDimension << std::endl;

  RIndexType rindex;
  rindex[0] = 0.5*( rsize[0] - m_TileDimension[1] );
  rindex[1] = 0.5*( rsize[1] - m_TileDimension[0] );
  rroi.SetIndex( rindex );

  //std::cout << rindex << std::endl;

  rsize[0] = m_TileDimension[1];
  rsize[1] = m_TileDimension[0];
  rroi.SetSize( rsize );

  //std::cout << rsize << std::endl;

  ROIFilterPointer roiFilter = ROIFilterType::New();
  roiFilter->SetRegionOfInterest( rroi );
  roiFilter->SetInput( currentImage );
  roiFilter->Update();
  RImagePointer tempImage = roiFilter->GetOutput();
  tempImage->DisconnectPipeline();

  RIteratorType It( tempImage, tempImage->GetLargestPossibleRegion() );
  It.GoToBegin();
  while( !It.IsAtEnd() )
  {
    double p = It.Get();
    if ( p < m_CorrectionThreshold )
    {
      It.Set( 0 );
    }
    else
    {
      It.Set( p - m_CorrectionThreshold );
    }
    ++It;
  }

  RescaleFilterPointer rescale = RescaleFilterType::New();
  rescale->SetInput( tempImage );
  rescale->SetOutputMinimum( 0.0 );
  rescale->SetOutputMaximum( 1.0 );
  rescale->Update();

  m_CorrectionImage = rescale->GetOutput();
  m_CorrectionImage->DisconnectPipeline();

  /*
  CastFilterPointer caster = CastFilterType::New();
  caster->SetInput( tempImage );
  caster->Update();

  WriterPointer writer = WriterType::New();
  writer->SetInput( caster->GetOutput() );
  writer->SetFileName( "/Users/kishoremosaliganti/CorrectionImage.tif" );
  writer->Update();
  */
}


/** Print Self information */
template<class TInputImage >
void
StitchingSharedData< TInputImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

} // end namespace

#endif
