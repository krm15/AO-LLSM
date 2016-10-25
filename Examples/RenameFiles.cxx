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

#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "itkDirectory.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"

int main ( int argc, char* argv[] )
{
  if ( argc < 2 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " iImageDir <rescaleFactors> <channelID>";
    std::cerr << "" << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int Dimension = 3;
  typedef itk::Directory DirectoryType;

  typedef double PixelType;
  typedef itk::Image< PixelType, Dimension > InputImageType;
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::ImageRegionIterator< InputImageType > IteratorType;
  typedef unsigned char OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter< InputImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  std::string m_TileDirectory = argv[1];

  bool m_RescaleImage = true;
  unsigned int histogramSize = 150000000;
  unsigned int discreteLimit = 255;//65535
  double percentileLimit = 0.001;
  std::vector< unsigned int > histogram;
  unsigned int p;
  double m_ScalingFactor = 1.0;
  unsigned int channelOfInterest = 0;

  if ( argc > 2 )
  {
    channelOfInterest = atoi( argv[2] );
  }

  if ( argc > 3 )
  {
    m_RescaleImage = false;
    m_ScalingFactor = atof( argv[3] );
  }

  DirectoryType::Pointer directory = DirectoryType::New();
  directory->Load( m_TileDirectory.c_str() );

  std::string searchString1 = "nm";
  std::string searchString2 = "stack";
  std::string searchString3 = "_ch";
  std::string searchString4 = "t_";
  for ( unsigned int m = 0; m < directory->GetNumberOfFiles(); m++)
  {
    //std::cout << "m: " << m << std::endl;
    std::string filename = directory->GetFile( m );
    if ( ( filename.find( searchString1 ) != std::string::npos ) &&
         ( filename.find( searchString2 ) != std::string::npos ) &&
         ( filename.find( searchString3 ) != std::string::npos ) &&
         ( filename.find( searchString4 ) != std::string::npos ) )
    {
      std::cout << filename << std::endl;

      unsigned int lastindex1 = 0; //filename.find_last_of( "." );
      unsigned int lastindex2 = filename.find_last_of( "." );
      //std::string path = filename.substr( 0, lastindex1+1 );
      std::string rawname = filename.substr( lastindex1+1, lastindex2 - lastindex1 -1 );

      //std::cout << "Path: " << path << std::endl;
      //std::cout << "Rawname: " << rawname << std::endl;

      unsigned int pos = rawname.find( searchString1 );
      std::string m_ChannelName = rawname.substr( pos-3, 5 );

      pos = rawname.find( searchString2 );
      std::string stack = rawname.substr( pos+5, 4 );
      //std::cout << "Stack: " << stack << std::endl;


      pos = rawname.find( searchString3 );
      std::string channel = rawname.substr( pos+3, 1 );
      //std::cout << "Channel: " << channel << std::endl;


      pos = rawname.find( searchString4 );
      std::string tp = rawname.substr( pos-4, 4 );
      //std::cout << "Time: " << tp << std::endl;


      if ( channelOfInterest == atoi( channel.c_str() ) )
      {
        std::string inputImageName = m_TileDirectory + filename;
        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName ( inputImageName.c_str() );
        reader->Update();
        InputImageType::Pointer inputImg = reader->GetOutput();

        if ( m_RescaleImage )
        {
          // Rescale the file
          histogram.resize( histogramSize, 0 );
          IteratorType It( inputImg, inputImg->GetLargestPossibleRegion() );
          It.GoToBegin();
          while( !It.IsAtEnd() )
          {
            p = static_cast<unsigned int>( It.Get() );
            if ( p >= histogramSize )
            {
              p = histogramSize-1;
            }
            histogram[p]++;
            ++It;
          }

          unsigned int totalPixelCount = inputImg->GetLargestPossibleRegion().GetNumberOfPixels();
          unsigned int topPercentOfPixels = percentileLimit * totalPixelCount;
          unsigned int maxIndex = histogramSize - 1;
          unsigned int cumsum = 0;
          while( ( maxIndex > 0 ) && ( cumsum < topPercentOfPixels ) )
          {
            cumsum += histogram[maxIndex];
            maxIndex--;
          }

          if ( maxIndex > 0 )
          {
            m_ScalingFactor = double(discreteLimit)/( (double)maxIndex );
          }
          std::cout << "Scaling Factor: " << m_ScalingFactor << std::endl;
        }

        if ( m_ScalingFactor < 0.999 )
        {
          IteratorType It( inputImg, inputImg->GetLargestPossibleRegion() );
          It.GoToBegin();
          while( !It.IsAtEnd() )
          {
            p = It.Get()*m_ScalingFactor;
            if ( p > discreteLimit )
            {
              p = discreteLimit;
            }

            It.Set( p );
            ++It;
          }
        }

        // Cast it to a different pixel type
        CastFilterType::Pointer caster = CastFilterType::New();
        caster->SetInput( inputImg );
        caster->Update();

        std::string oFilename = m_TileDirectory + m_ChannelName
                + "_" + stack + "stack_" + channel + "ch_" + tp + "t.mha";
        std::cout << "Output filename: " << oFilename << std::endl;

        // Write out the file
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName( oFilename.c_str() );
        writer->SetInput( caster->GetOutput() );
        writer->UseCompressionOn();
        writer->Update();
      }
    }
  }

  return EXIT_SUCCESS;
  }
