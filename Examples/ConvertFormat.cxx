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
#include "itkImageRegionIteratorWithIndex.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkSettingsInfoExtractionFilter.h"

int main ( int argc, char* argv[] )
{
  if ( argc < 5 )
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " iInputSettingsDir iInputImageDir oOutputImageDir ";
    std::cerr << "iChannelNumber iTimePoint" << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int Dimension = 3;
  typedef std::vector< std::string > StringVectorType;
  typedef std::vector< double > DoubleVectorType;
  typedef vnl_matrix< double > vnlMatrixType;
  typedef vnl_vector< double > vnlVectorType;
  typedef itk::Directory DirectoryType;
  typedef std::vector< std:: vector< std::vector< std::string > > > StringArray3DType;

  typedef unsigned short PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  typedef itk::PermuteAxesImageFilter< ImageType > PermuteAxesFilterType;

  typedef itk::SettingsInfoExtractionFilter< double, ImageType > SettingsFilterType;

  typedef unsigned short OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter< ImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  typedef ImageType::SpacingType SpacingType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::RegionType RegionType;
  typedef ImageType::PointType PointType;
  typedef SizeType::SizeValueType SizeValueType;

  SettingsFilterType::Pointer settingsReader = SettingsFilterType::New();
  settingsReader->SetSettingsDirectory( argv[1] );
  settingsReader->SetTileDirectory( argv[2] );
  settingsReader->SetChannelNumber( atoi(argv[4]) );
  settingsReader->SetTimePoint( atoi(argv[5]) );
  settingsReader->Read();

  StringVectorType m_SettingName = settingsReader->GetSettingFieldName();
  DoubleVectorType m_SettingValue = settingsReader->GetSettingFieldValue(  );

  // Setup the dimensions of the largest stitched image
  unsigned int numOfTiles = settingsReader->GetNumberOfTiles();
  unsigned int *tileNumber;
  tileNumber = settingsReader->GetTileNumber();

  double *tileSize;
  tileSize = settingsReader->GetTileSize();

  SizeType tilePixelDimension = settingsReader->GetTileDimension();
  SpacingType spacing = settingsReader->GetTileSpacing();

  std::cout << "Number of tiles " << numOfTiles << std::endl;
  std::cout << "Tile number" << std::endl;
  std::cout << tileNumber[0] << ' ' << tileNumber[1] << ' '
                             << tileNumber[2] << std::endl;
  std::cout << " Tile size (um)" << std::endl;
  std::cout << tileSize[0] << ' ' << tileSize[1] << ' '
                           << tileSize[2] << std::endl;
  std::cout << "Tile pixel dimension" << std::endl;
  std::cout << tilePixelDimension << std::endl;
  std::cout << "Tile spacing" << std::endl;
  std::cout << spacing << std::endl;


  settingsReader->CreateStitchedImage();

  std::cout << std::endl;
  std::cout << "Stitched image origin" << std::endl;
  std::cout << settingsReader->GetStitchOrigin() << std::endl;
  std::cout << "Stitched image dimensions" << std::endl;
  std::cout << settingsReader->GetStitchDimension() << std::endl;
  std::cout << "Stitched image size" << std::endl;
  std::cout << settingsReader->GetStitchSize()[0] << ' ';
  std::cout << settingsReader->GetStitchSize()[1] << ' ';
  std::cout << settingsReader->GetStitchSize()[2] << std::endl;


  // Read all the files in the input directory of type ch and at timePoint
  std::string filename;
  itk::FixedArray<unsigned int, 3> axesOrder;
  axesOrder[0] = 1;
  axesOrder[1] = 0;
  axesOrder[2] = 2;

  for( unsigned int i = 0; i < tileNumber[0]; i++ )
  {
    for( unsigned int j = 0; j < tileNumber[1]; j++ )
    {
      for( unsigned int k = 0; k < tileNumber[2]; k++ )
      {
        filename = settingsReader->GetTileFileNameArray()[i][j][k];
        std::cout << "Filename: " << filename << std::endl;
        std::stringstream  filename3;

        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName ( filename.c_str() );
        reader->Update();

//        PermuteAxesFilterType::Pointer pAFilter = PermuteAxesFilterType::New();
//        pAFilter->SetInput( reader->GetOutput() );
//        pAFilter->SetOrder( axesOrder );
//        pAFilter->Update();
//        ImageType::Pointer pImage = pAFilter->GetOutput();

//        ImageType::Pointer currentImage = ImageType::New();
//        currentImage->SetOrigin( pImage->GetOrigin() );
//        currentImage->SetSpacing( pImage->GetSpacing() );
//        currentImage->SetRegions( pImage->GetLargestPossibleRegion() );
//        currentImage->Allocate();

//        IteratorType pIt( pImage, pImage->GetLargestPossibleRegion() );
//        IteratorType cIt( currentImage, currentImage->GetLargestPossibleRegion() );
//        while(!pIt.IsAtEnd())
//        {
//          cIt.Set( pIt.Get() );
//          ++pIt;
//          ++cIt;
//        }

//        CastFilterType::Pointer caster = CastFilterType::New();
//        caster->SetInput( currentImage );//rescale->GetOutput()
//        caster->Update();

        unsigned int lastindex1 = filename.find_last_of( "/" );
        unsigned int lastindex2 = filename.find_last_of( "." );
        std::string rawname = filename.substr( lastindex1+1, lastindex2 - lastindex1 -1 );
        std::cout << rawname << std::endl;
        filename3 << argv[3] << rawname << ".mha";
        std::cout << filename3.str().c_str() << std::endl;

        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName( filename3.str().c_str() );
        writer->SetInput( reader->GetOutput() );
        //writer->UseCompressionOn();
        writer->Update();
      }
    }
  }

  return EXIT_SUCCESS;
  }
