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

#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>

#include "itkDirectory.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"

int main ( int argc, char* argv[] )
{
  const unsigned int Dimension = 2;
  typedef double PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::Directory DirectoryType;
  typedef itk::ImageRegionIterator< ImageType > IteratorType;

  typedef ImageType::SpacingType SpacingType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::RegionType RegionType;
  typedef ImageType::PointType PointType;
  typedef SizeType::SizeValueType SizeValueType;

  std::string searchCH = "_ch";

  DirectoryType::Pointer directory = DirectoryType::New();
  directory->Load( argv[1] );

  std::string filename;
  std::stringstream searchString;
  searchString << "_MIP_z.tif";

  unsigned int histogramSize = 500000;
  std::vector< unsigned int > histogram;
  unsigned int m_TrueCountOfTiles = 0;
  unsigned int totalPixelCount = 0;
  PixelType p;
  for ( unsigned int m = 0; m < directory->GetNumberOfFiles(); m++)
  {
    //std::cout << "m: " << m << std::endl;
    filename = directory->GetFile( m );
    if ( filename.find( searchString.str() ) != std::string::npos )
    {
      std::ifstream infile( filename.c_str() );
      if ( infile )
      {
       infile.close();

       std::cout << filename << std::endl;
       ReaderType::Pointer reader = ReaderType::New();
       reader->SetFileName ( filename.c_str() );

       try
       {
         reader->Update();
       }
       catch( itk::ExceptionObject & err )
       {
         std::cerr << "ExceptionObject caught !" << std::endl;
         std::cerr << err << std::endl;
       }

       ImageType::Pointer img = reader->GetOutput();
       totalPixelCount += img->GetLargestPossibleRegion().GetNumberOfPixels();

       IteratorType It( img, img->GetLargestPossibleRegion() );
       It.GoToBegin();
       while( !It.IsAtEnd() )
       {
         p = static_cast<unsigned int>(It.Get());
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

  double rescaleFactor = 1.0;
  if ( totalPixelCount > 0 )
  {
    unsigned int topPercentOfPixels = 0.03 * totalPixelCount;
    unsigned int max = histogramSize - 1;
    unsigned int cumsum = 0;
    while( ( max > 0 ) && ( cumsum < topPercentOfPixels ) )
    {
      cumsum += histogram[max];
      max--;
    }

    if ( max > 0 )
    {
     rescaleFactor = double(65535)/( (double)max );
    }
  }
  std::cout << "Scaling Factor: " << rescaleFactor << std::endl;

  std::ofstream outputFile( argv[2] );
  outputFile << rescaleFactor << std::endl;
  outputFile.close();

  return EXIT_SUCCESS;
}
