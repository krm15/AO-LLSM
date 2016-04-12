/*=========================================================================
  Author: $Author: arnaudgelas $  // Author of last commit
  Version: $Rev: 567 $  // Revision of last commit
  Date: $Date: 2009-08-17 11:47:32 -0400 (Mon, 17 Aug 2009) $  // Date of last commit
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
/*=========================================================================
  Author: $Author: arnaudgelas $  // Author of last commit
  Version: $Rev: 567 $  // Revision of last commit
  Date: $Date: 2009-08-17 11:47:32 -0400 (Mon, 17 Aug 2009) $  // Date of last commit
=========================================================================*/

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRichardsonLucyDeconvolutionImageFilter.h"

int main(int argc, char* argv[])
{
  if ( argc < 5 )
    {
    std::cerr << "Usage: " << argv[0]
              << " iInputImage iKernelImage oOutputImage iterations"
              << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef float                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  typedef itk::RichardsonLucyDeconvolutionImageFilter< ImageType > DeconvolutionFilterType;

  unsigned int iterations = static_cast< unsigned int >( atoi( argv[4] ) );

  ReaderType::Pointer inputReader = ReaderType::New();
  inputReader->SetFileName( argv[1] );
  inputReader->Update();

  ReaderType::Pointer kernelReader = ReaderType::New();
  kernelReader->SetFileName( argv[2] );
  kernelReader->Update();

  DeconvolutionFilterType::Pointer deconvolutionFilter = DeconvolutionFilterType::New();
  deconvolutionFilter->SetInput( inputReader->GetOutput() );
  deconvolutionFilter->SetKernelImage( kernelReader->GetOutput() );
  deconvolutionFilter->NormalizeOn();
  deconvolutionFilter->SetNumberOfIterations( iterations );

  // Write the deconvolution result
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( deconvolutionFilter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
