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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"

int main ( int argc, char* argv[] )
  {
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " input output <startX> <startY> <startZ> <endX> <endY> <endZ>" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef itk::Image< unsigned char, Dimension > InputImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::ImageFileWriter< InputImageType > WriterType;
  typedef itk::ImageRegionIterator< InputImageType > IteratorType;
  typedef itk::RegionOfInterestImageFilter< InputImageType,
    InputImageType > ROIFilterType;
  typedef InputImageType::RegionType RegionType;
  typedef InputImageType::IndexType IndexType;
  typedef IndexType::IndexValueType IndexValueType;
  typedef InputImageType::SizeType SizeType;
  typedef SizeType::SizeValueType SizeValueType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();

  RegionType region = reader->GetOutput()->GetLargestPossibleRegion();
  SizeType size = region.GetSize();

  IndexType start, end;
  start = region.GetIndex();

  if (argc > 3)
    {
    start[0] = static_cast<IndexValueType>( atoi( argv[3] ) );
    start[1] = static_cast<IndexValueType>( atoi( argv[4] ) );
    start[2] = static_cast<IndexValueType>( atoi( argv[5] ) );

    end[0] = static_cast<IndexValueType>( atoi( argv[6] ) );
    end[1] = static_cast<IndexValueType>( atoi( argv[7] ) );
    end[2] = static_cast<IndexValueType>( atoi( argv[8] ) );

    for(unsigned int i=0; i<Dimension; i++)
      {
      size[i] = end[i]-start[i]+1;
      }
    }
  std::cout << size << std::endl;

  InputImageType::RegionType subRegion;
  subRegion.SetIndex( start );
  subRegion.SetSize( size );

  ROIFilterType::Pointer roi = ROIFilterType::New();
  roi->SetInput( reader->GetOutput() );
  roi->SetRegionOfInterest( subRegion );
  roi->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput ( roi->GetOutput() );
  writer->SetFileName ( argv[2] );

  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject e )
    {
    std::cerr << "Error: " << e << std::endl;
    }

  return EXIT_SUCCESS;
  }
