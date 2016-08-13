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
#include "itkRichardsonLucyDeconvolutionImageFilter.h"
#include "itkStitchingSharedData.h"
#include "anyoption.h"

int main ( int argc, char* argv[] )
{
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
  typedef itk::RichardsonLucyDeconvolutionImageFilter< ImageType > DeconvolutionFilterType;

  typedef unsigned short OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter< ImageType, OutputImageType > CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  typedef itk::StitchingSharedData< ImageType > SharedDataType;

  typedef ImageType::SpacingType SpacingType;
  typedef ImageType::SizeType SizeType;
  typedef ImageType::IndexType IndexType;
  typedef ImageType::RegionType RegionType;
  typedef ImageType::PointType PointType;
  typedef SizeType::SizeValueType SizeValueType;


  /* 1. CREATE AN OBJECT */
  AnyOption *opt = new AnyOption();

  /* 2. SET PREFERENCES  */
  //opt->noPOSIX(); /* do not check for POSIX style character options */
  //opt->setVerbose(); /* print warnings about unknown options */
  //opt->autoUsagePrint(true); /* print usage for bad options */

  /* 3. SET THE USAGE/HELP   */
  opt->addUsage( "" );
  opt->addUsage( "Usage: " );
  opt->addUsage( "" );
  opt->addUsage( " iSettings file directory " );
  opt->addUsage( " iTile directories " );
  opt->addUsage( " Output directory " );
  opt->addUsage( " -h   --help    Prints this help " );
  opt->addUsage( " -c   --channel 0   (default) channel value" );
  opt->addUsage( " -t   --time    0   (default) timepoint" );
  opt->addUsage( " -d   --deconv  ~/  (default) deconvolve tiles based on PSF" );
  opt->addUsage( "" );

  /* 4. SET THE OPTION STRINGS/CHARACTERS */

  /* by default all  options will be checked on the command line
    and from option/resource file */

  /* a flag (takes no argument), supporting long and short form */
  opt->setFlag(  "help",  'h' );

  /* an option (takes an argument), supporting long and short form */
  opt->setOption(  "channel", 'c' );
  opt->setOption(  "time",    't' );
  opt->setOption(  "deconv",  'd' );

  /* 5. PROCESS THE COMMANDLINE AND RESOURCE FILE */
  /* read options from a  option/resource file with ':'
  separated options or flags, one per line */

  opt->processFile( ".options" );
  opt->processCommandArgs( argc, argv );

  if( ! opt->hasOptions())
  {
    opt->printUsage();
    delete opt;
    return EXIT_FAILURE;
  }

  /* 6. GET THE VALUES */
  if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) || ( opt->getArgc() < 3 ) )
  {
    opt->printUsage();
    delete opt;
    return EXIT_FAILURE;
  }

  unsigned int ch = 0;
  unsigned int tp = 0;
  unsigned int numOfIterations = 15;
  bool deconv = false;
  std::string PSFImagePath = "~/";

  if( opt->getValue( 'c' ) != NULL  || opt->getValue( "channel" ) != NULL  )
  {
    ch = atoi( opt->getValue( 'c' ) );
  }
  if( opt->getValue( 't' ) != NULL  || opt->getValue( "time" ) != NULL  )
  {
    tp = atoi( opt->getValue( 't' ) );
  }
  if( opt->getValue( 'd' ) != NULL  || opt->getValue( "deconv" ) != NULL  )
  {
    PSFImagePath = opt->getValue( 'd' );
    deconv = true;
  }

  SharedDataType::Pointer m_SharedData = SharedDataType::New();

  SettingsFilterType::Pointer settingsReader = SettingsFilterType::New();
  settingsReader->SetSettingsDirectory( argv[1] );
  settingsReader->SetTileDirectory( argv[2] );
  settingsReader->SetChannelNumber( ch );
  settingsReader->SetTimePoint( tp );
  settingsReader->SetSharedData( m_SharedData );
  settingsReader->Read();

  // Setup the dimensions of the largest stitched image
  unsigned int numOfTiles = settingsReader->GetNumberOfTiles();
  IndexType tileNumber = m_SharedData->m_TileNumber;
  PointType tileSize = m_SharedData->m_TileSize;
  SizeType tilePixelDimension = m_SharedData->m_TileDimension;
  SpacingType spacing = m_SharedData->m_TileSpacing;
  PointType tileOverlap = m_SharedData->m_TileOverlap;

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

  ImageType::Pointer kernelImage;
  if ( deconv )
  {
    ReaderType::Pointer kernelReader = ReaderType::New();
    kernelReader->SetFileName( PSFImagePath );
    kernelReader->Update();
    kernelImage = kernelReader->GetOutput();
    kernelImage->DisconnectPipeline();
  }


  // Read all the files in the input directory of type ch and at timePoint
  std::string filename;
  itk::FixedArray<unsigned int, 3> axesOrder;
  axesOrder[0] = 1;
  axesOrder[1] = 0;
  axesOrder[2] = 2;

  for( unsigned int k = 0; k < tileNumber[2]; k++ )
  {
    for( unsigned int i = 0; i < tileNumber[0]; i++ )
    {
       for( unsigned int j = 0; j < tileNumber[1]; j++ )
       {
         filename = settingsReader->GetSharedData()->GetTileFileNameArray()[i][j][k];

         //std::cout << "Filename: " << filename << std::endl;
         if ( !filename.empty() )
         {
           ReaderType::Pointer reader = ReaderType::New();
           reader->SetFileName ( filename.c_str() );
           reader->Update();

//        PermuteAxesFilterType::Pointer pAFilter = PermuteAxesFilterType::New();
//        pAFilter->SetInput( reader->GetOutput() );
//        pAFilter->SetOrder( axesOrder );
//        pAFilter->Update();
           ImageType::Pointer pImage = reader->GetOutput();

           ImageType::Pointer currentImage = ImageType::New();
           currentImage->SetOrigin( pImage->GetOrigin() );
           currentImage->SetSpacing( pImage->GetSpacing() );
           currentImage->SetRegions( pImage->GetLargestPossibleRegion() );
           currentImage->Allocate();

           IteratorType pIt( pImage, pImage->GetLargestPossibleRegion() );
           IteratorType cIt( currentImage, currentImage->GetLargestPossibleRegion() );
           while(!pIt.IsAtEnd())
           {
             cIt.Set( pIt.Get() );
             ++pIt;
             ++cIt;
           }

           ImageType::Pointer deconvImage;
           if ( deconv )
           {
             DeconvolutionFilterType::Pointer deconvolutionFilter = DeconvolutionFilterType::New();
             deconvolutionFilter->SetInput( currentImage );
             deconvolutionFilter->SetKernelImage( kernelImage );
             deconvolutionFilter->NormalizeOn();
             deconvolutionFilter->SetNumberOfIterations( numOfIterations );
             deconvolutionFilter->Update();
             deconvImage = deconvolutionFilter->GetOutput();
             deconvImage->DisconnectPipeline();
           }
           else
           {
             deconvImage = currentImage;
           }

           CastFilterType::Pointer caster = CastFilterType::New();
           caster->SetInput( deconvImage );//rescale->GetOutput()
           caster->Update();

           std::stringstream  filename3;
           unsigned int lastindex1 = filename.find_last_of( "/" );
           unsigned int lastindex2 = filename.find_last_of( "." );
           std::string rawname = filename.substr( lastindex1+1, lastindex2 - lastindex1 -1 );
           std::cout << rawname << std::endl;
           filename3 << argv[3] << rawname << ".mha";
           std::cout << filename3.str().c_str() << std::endl;

           WriterType::Pointer writer = WriterType::New();
           writer->SetFileName( filename3.str().c_str() );
           writer->SetInput( caster->GetOutput() );
           writer->Update();
        }
      }
    }
  }

  return EXIT_SUCCESS;
  }
