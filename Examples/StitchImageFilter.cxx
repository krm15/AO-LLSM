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
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericSeriesFileNames.h"

int main ( int argc, char* argv[] )
  {
  if ( argc < 3 )
    {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " iInputImageFormat oOutputImageFormat startX startY startZ endX endY endZ sizeX sizeY sizeZ" << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef unsigned short PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  typedef itk::NumericSeriesFileNames NameGeneratorType;

  unsigned int startX = atoi(argv[3]);
  unsigned int startY = atoi(argv[4]);
  unsigned int startZ = atoi(argv[5]);
  unsigned int stopX = atoi(argv[6]);
  unsigned int stopY = atoi(argv[7]);
  unsigned int stopZ = atoi(argv[8]);

  // Allocated output image
  ImageType::PointType origin;
  ImageType::SizeType size;
  size[0] = atoi(argv[9]);
  size[1] = atoi(argv[10]);
  size[2] = atoi(argv[11]);

  ImageType::IndexType index;
  for(unsigned int i = 0; i < Dimension; i++)
  {
      origin[i] = 0.0;
      index[i] = 0;
  }

  ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( index );

  ImageType::Pointer outputImage = ImageType::New();
  outputImage->SetOrigin( origin );
  outputImage->SetRegions( region );
  outputImage->Allocate();
  outputImage->FillBuffer( 0 );

  for(unsigned int i = 0;; i < numberOfFiles; i++ )
  {
    //Read each image
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName ( nameGeneratorOutput->GetFileNames()[i] );
    reader->Update();
    ImageType::Pointer currentImage = reader->GetOutput();
    currentImage->DisconnectPipeline();

    //Fill output image

    subRegion

    IteratorType cIt( currentImage, region );
    IteratorType oIt( output, subRegion );

    cIt.GoToBegin();
    oIt.GoToBegin();
    while( !cIt.IsAtEnd() )
    {
      oIt.Set( cIt.Get() );
      ++cIt;
    }
  }

  // Write output image
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput ( outputImage );
  writer->SetFileName ( argv[6] );

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
