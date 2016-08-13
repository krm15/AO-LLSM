/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 1550 $  // Revision of last commit
  Date: $Date: 2010-06-06 23:50:34 -0400 (Sun, 06 Jun 2010) $  // Date of last commit
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

#ifndef __itkExtractStitchTransformImageFilter_h
#define __itkExtractStitchTransformImageFilter_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <vector>

#include <cstring>
#include <strstream>

#include <iostream>
#include <istream>
#include <streambuf>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <vnl_vector.h>
#include <vnl_matrix.h>

#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkStitchingSharedData.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkMultiThreader.h"

namespace itk
{
template < class TInputImage >
class ITK_EXPORT ExtractStitchTransformImageFilter : public Object
{
  public:
  typedef ExtractStitchTransformImageFilter Self;
  typedef Object                            Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Method for creation through object factory */
  itkNewMacro ( Self );

  /** Run-time type information */
  itkTypeMacro ( ExtractStitchTransformImageFilter, Object );

  typedef std::vector< std::string > StringVectorType;
  typedef std::vector< double > DoubleVectorType;

  typedef TInputImage ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef ImageFileReader< ImageType > ReaderType;
  typedef typename ReaderType::Pointer ReaderPointer;
  typedef itk::ImageRegionIterator< ImageType > IteratorType;

  typedef StitchingSharedData< ImageType > SharedDataType;
  typedef typename SharedDataType::Pointer SharedDataPointer;

  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::RegionType RegionType;
  typedef typename ImageType::PointType PointType;
  typedef typename SizeType::SizeValueType SizeValueType;
  typedef std::vector< IndexType > IndexVectorType;

  typedef RegionOfInterestImageFilter< ImageType, ImageType > ROIFilter3DType;
  typedef typename ROIFilter3DType::Pointer ROIFilter3DPointer;

  typedef MultiThreader ThreaderType;
  typedef typename ThreaderType::Pointer ThreaderPointer;

  void SetStepLength( DoubleVectorType& stepLength )
  {
    m_StepLength = stepLength;
  }

  void SetSearchRadius( DoubleVectorType& searchRadius )
  {
    m_SearchRadius = searchRadius;
  }

  itkGetObjectMacro( SharedData, SharedDataType );

  itkSetMacro( TileDirectory,       std::string );
  itkSetMacro( OffsetFilePath,      std::string );
  itkSetMacro( ZTileStart,          unsigned int );
  itkSetMacro( ZTileEnd,            unsigned int );
  itkSetMacro( NumberOfThreads,     unsigned int );
  itkSetObjectMacro( SharedData,    SharedDataType );

  void BeginRegister();

  protected:
  ExtractStitchTransformImageFilter();
  ~ExtractStitchTransformImageFilter(){}
  void OverlapRegion( ImagePointer A, ImagePointer B,
                      RegionType& rA, RegionType& rB );
  void ReadOffsetFile();
  void WriteOffsetFile();
  void RegisterTiles( float searchRadius, float stepLength, unsigned int ztile );
  static ITK_THREAD_RETURN_TYPE ThreaderCallback(void * arg);

  struct ThreadStruct
  {
    Self*                 Filter;
  };

  std::string   m_TileDirectory;
  std::string   m_OffsetFilePath;

  unsigned int      m_Dimension;
  unsigned int      m_ZTileStart;
  unsigned int      m_ZTileEnd;
  unsigned int      m_NumberOfThreads;
  SharedDataPointer m_SharedData;
  DoubleVectorType  m_SearchRadius;
  DoubleVectorType  m_StepLength;

  private:
    ExtractStitchTransformImageFilter ( Self& );   // intentionally not implemented
    void operator= ( const Self& );   // intentionally not implemented
  };

} /* namespace itk */

#include "itkExtractStitchTransformImageFilter.txx"
#endif
