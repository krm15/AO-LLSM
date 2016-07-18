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
#include "itkRichardsonLucyDeconvolutionImageFilter.h"

int main(int argc, char* argv[])
{
  if ( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0]
              << " iInputImage iKernelImage oOutputImage <iterations>"
              << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 3;
  typedef float                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType >  WriterType;
  typedef itk::RichardsonLucyDeconvolutionImageFilter< ImageType > DeconvolutionFilterType;

  unsigned int m_NumberOfIterations = 15;
  if ( argc > 4)
  {
    m_NumberOfIterations = static_cast< unsigned int >( atoi( argv[4] ) );
    std::cout << "Number of iterations: " << m_NumberOfIterations << std::endl;
  }

  ReaderType::Pointer inputReader = ReaderType::New();
  inputReader->SetFileName( argv[1] );
  inputReader->Update();
  ImageType::Pointer input = inputReader->GetOutput();
  input->DisconnectPipeline();
  std::cout << "Input image read" << std::endl;

  ImageType::SpacingType spacing;
  spacing[0] = spacing[1] = 0.09;
  spacing[2] = 0.2;

  input->SetSpacing( spacing );

  ReaderType::Pointer kernelReader = ReaderType::New();
  kernelReader->SetFileName( argv[2] );
  kernelReader->Update();
  ImageType::Pointer kernel = kernelReader->GetOutput();
  kernel->DisconnectPipeline();
  std::cout << "Kernel image read" << std::endl;

  spacing[0] = spacing[1] = 0.09;
  spacing[2] = 0.1;
  kernel->SetSpacing( spacing );

  std::cout << "Beginning deconvolution" << std::endl;
  DeconvolutionFilterType::Pointer deconvolutionFilter = DeconvolutionFilterType::New();
  deconvolutionFilter->SetInput( inputReader->GetOutput() );
  deconvolutionFilter->SetKernelImage( kernelReader->GetOutput() );
  deconvolutionFilter->NormalizeOn();
  deconvolutionFilter->SetNumberOfIterations( m_NumberOfIterations );
  deconvolutionFilter->SetNumberOfThreads( 2 );
  deconvolutionFilter->Update();
  std::cout << "Ending deconvolution" << std::endl;

  // Write the deconvolution result
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[3] );
  writer->SetInput( deconvolutionFilter->GetOutput() );
  writer->Update();

  return EXIT_SUCCESS;
}
