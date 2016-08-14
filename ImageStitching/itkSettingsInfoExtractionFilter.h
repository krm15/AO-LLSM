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

#ifndef __itkSettingsInfoExtractionFilter_h
#define __itkSettingsInfoExtractionFilter_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include <vector>

#include <cstring>
#include <iostream>
#include <istream>
#include <streambuf>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <vnl_vector.h>
#include <vnl_matrix.h>
#include "itkDirectory.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkStitchingSharedData.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

#include "itkFixedArray.h"
#include "itk_tiff.h"

#include "itkImageFileWriter.h"


namespace itk
{
template < class TValueType, class TInputImage >
class ITK_EXPORT SettingsInfoExtractionFilter : public Object
{
  public:
  typedef SettingsInfoExtractionFilter Self;
  typedef Object                  Superclass;
  typedef SmartPointer< Self >         Pointer;
  typedef SmartPointer< const Self >   ConstPointer;

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Method for creation through object factory */
  itkNewMacro ( Self );

  /** Run-time type information */
  itkTypeMacro ( SettingsInfoExtractionFilter, Object );

  typedef TValueType ValueType;

  typedef std::vector< std::string > StringVectorType;
  typedef std::vector< ValueType > DoubleVectorType;
  typedef vnl_matrix< ValueType > vnlMatrixType;
  typedef vnl_vector< ValueType > vnlVectorType;
  typedef Directory DirectoryType;
  typedef typename DirectoryType::Pointer DirectoryPointer;

  typedef TInputImage ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef ImageFileReader< ImageType > ReaderType;
  typedef typename ReaderType::Pointer ReaderPointer;

  typedef Image<double, 2> RImageType;
  typedef typename RImageType::Pointer RImagePointer;

  typedef ImageFileReader< RImageType > RReaderType;
  typedef typename RReaderType::Pointer RReaderPointer;
  typedef ImageRegionIterator< RImageType > RIteratorType;

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

  typedef MetaDataDictionary DictionaryType;
  typedef MetaDataObject< std::string > MetaDataStringType;

  void Read();
  void UpdateFileNameLookup( std::istream& os );
  void CreateStitchedImage();
  void AllocateROI();

  itkGetObjectMacro( StitchedImage, ImageType );
  itkGetObjectMacro( SharedData, SharedDataType );

  itkGetConstMacro( SettingFieldName,   StringVectorType );
  itkGetConstMacro( SettingFieldValue,  StringVectorType );
  itkGetConstMacro( NumberOfTiles,      unsigned int );
  itkGetConstMacro( MinimumStart,       PointType );
  itkGetConstMacro( MaximumEnd,         PointType );
  itkGetConstMacro( StitchOrigin,       PointType );
  itkGetConstMacro( StitchDimension,    SizeType );
  itkGetConstMacro( StitchIndex,        IndexType );
  itkGetConstMacro( StitchSize,         PointType );
  itkGetConstMacro( StitchRegion,       RegionType );
  itkGetConstMacro( ChannelName,        std::string );
  itkGetConstMacro( ScopeName,          std::string );
  itkGetConstMacro( SampleScan,         bool );

  itkSetMacro( SettingsDirectory,   std::string );
  itkSetMacro( TileDirectory,       std::string );
  itkSetMacro( OffsetFilePath,      std::string );
  itkSetMacro( ChannelPrefix,       std::string );
  itkSetMacro( ChannelNumber,       unsigned int );
  itkSetMacro( TimePoint,           unsigned int );
  itkSetObjectMacro( SharedData, SharedDataType );

  protected:
  SettingsInfoExtractionFilter();
  ~SettingsInfoExtractionFilter(){}
  void UpdateTileCoverage( std::istream& os );
  void TransformCoordinateAxes();
  void ReadTileInfo( std::istream& os );
  void ReadOffsetFile();
  std::istream& safeGetline(std::istream& is, std::string& t);

  std::string   m_TileDirectory;
  std::string   m_OffsetFilePath;
  std::string   m_ChannelName;
  std::string   m_SampleName;
  unsigned int  m_ChannelNumber;
  unsigned int  m_TimePoint;
  std::string   m_ChannelPrefix;

  unsigned int      m_Dimension;
  std::string       m_SettingsDirectory;
  StringVectorType  m_SettingFieldName;
  StringVectorType  m_SettingFieldValue;
  bool              m_SampleScan;
  std::string       m_ScopeName;

  unsigned int      m_NumberOfTiles;

  PointType m_MinimumStart;
  PointType m_MaximumEnd;

  vnlMatrixType m_TileInfoValue;
  vnlMatrixType m_TransformedTileInfoValue;

  ImagePointer  m_StitchedImage;
  PointType     m_StitchSize;
  PointType     m_StitchOrigin;
  SizeType      m_StitchDimension;
  IndexType     m_StitchIndex;
  RegionType    m_StitchRegion;

  FixedArray<unsigned int, ImageDimension> tileAxesOrder;

  SharedDataPointer m_SharedData;

  private:
    SettingsInfoExtractionFilter ( Self& );   // intentionally not implemented
    void operator= ( const Self& );   // intentionally not implemented
  };

} /* namespace itk */

#include "itkSettingsInfoExtractionFilter.txx"
#endif
