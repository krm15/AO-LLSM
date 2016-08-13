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
#ifndef __itkFillROIImageFilter_h
#define __itkFillROIImageFilter_h

#include <fstream>
#include "itkInPlaceImageFilter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkRichardsonLucyDeconvolutionImageFilter.h"
#include "itkStitchingSharedData.h"

namespace itk
{

template< class TInputImage >
class ITK_EXPORT FillROIImageFilter : public InPlaceImageFilter<TInputImage>
{
public:
  typedef FillROIImageFilter                Self;
  typedef InPlaceImageFilter<TInputImage>   Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( FillROIImageFilter, InPlaceImageFilter );

  itkStaticConstMacro ( ImageDimension, unsigned int, TInputImage::ImageDimension );

  /** Input Image typedef */
  typedef TInputImage                        ImageType;
  typedef typename ImageType::Pointer        ImagePointer;
  typedef typename ImageType::ConstPointer   ImageConstPointer;
  typedef typename ImageType::IndexType      IndexType;
  typedef typename IndexType::IndexValueType IndexValueType;
  typedef typename ImageType::PixelType      PixelType;
  typedef typename ImageType::SizeType       SizeType;
  typedef typename SizeType::SizeValueType   SizeValueType;
  typedef typename ImageType::RegionType     RegionType;
  typedef typename ImageType::SpacingType    SpacingType;
  typedef typename ImageType::PointType      PointType;
  typedef typename PointType::CoordRepType   CoordType;

  typedef ImageRegionIterator< ImageType > IteratorType;
  typedef ImageRegionIteratorWithIndex< ImageType > IndexIteratorType;

  typedef Image<unsigned short, ImageDimension> DoubleImageType; // double
  typedef ImageFileReader< DoubleImageType >    ReaderType;
  typedef typename ReaderType::Pointer          ReaderPointer;

  typedef CastImageFilter< DoubleImageType, ImageType > CastFilterType;
  typedef typename CastFilterType::Pointer CastFilterPointer;

  typedef PermuteAxesImageFilter< ImageType > PermuteAxesFilterType;
  typedef typename PermuteAxesFilterType::Pointer PermuteAxesFilterPointer;

  typedef RichardsonLucyDeconvolutionImageFilter< ImageType > DeconvolutionFilterType;
  typedef typename DeconvolutionFilterType::Pointer DeconvolutionFilterPointer;

  typedef RegionOfInterestImageFilter< ImageType, ImageType > ROIFilter3DType;
  typedef typename ROIFilter3DType::Pointer ROIFilter3DPointer;
  
  typedef ImageFileWriter< ImageType > WriterType;
  typedef typename WriterType::Pointer WriterPointer;

  typedef Image< double, 2 > RImageType;
  typedef typename RImageType::Pointer RImagePointer;
  typedef ImageFileReader< RImageType > RReaderType;
  typedef typename RReaderType::Pointer RReaderPointer;
  typedef ImageRegionIterator< RImageType > RIteratorType;

  typedef double ValueType;
  typedef std::vector< std::string > StringVectorType;
  typedef std::vector< std:: vector< std::vector< std::string > > > StringArray3DType;
  typedef StitchingSharedData< ImageType > SharedDataType;
  typedef typename SharedDataType::Pointer SharedDataPointer;

  itkSetObjectMacro( SharedData, SharedDataType );

  void SetZTile( unsigned int& ztile  )
  {
    m_ZTile = ztile;
    m_SingleZFill = true;
  }

protected:
  FillROIImageFilter();
  ~FillROIImageFilter(){}
  void OverlapRegion( ImagePointer A, ImagePointer B,
                      RegionType& rA, RegionType& rB );
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void BeforeThreadedGenerateData();
  virtual void AfterThreadedGenerateData(){}
  virtual void ThreadedGenerateData(const RegionType & windowRegion,
                                    ThreadIdType threadId);

  ImagePointer ExtractCorrectedAndFlippedTile( std::string& filename );

  SharedDataPointer m_SharedData;
  IndexType m_ScanStart;
  IndexType m_ScanEnd;
  unsigned int m_NumOfValidThreads;
  unsigned int      m_ZTile;
  bool m_SingleZFill;

private:
  FillROIImageFilter(const Self&) {}
  void operator=(const Self&) {}
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFillROIImageFilter.txx"
#endif

#endif
