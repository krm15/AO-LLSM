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

#include <iostream>
#include <vector>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vnl_vector.h>
#include <vnl_matrix.h>
#include "itkDirectory.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkPermuteAxesImageFilter.h"

namespace itk
{
template < class TValueType, class TInputImage >
class ITK_EXPORT SettingsInfoExtractionFilter : public LightObject
{
  public:
  typedef SettingsInfoExtractionFilter Self;
  typedef LightObject                  Superclass;
  typedef SmartPointer< Self >         Pointer;
  typedef SmartPointer< const Self >   ConstPointer;

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Method for creation through object factory */
  itkNewMacro ( Self );

  /** Run-time type information */
  itkTypeMacro ( SettingsInfoExtractionFilter, LightObject );

  typedef TValueType ValueType;

  typedef std::vector< std::string > StringVectorType;
  typedef std::vector< ValueType > DoubleVectorType;
  typedef vnl_matrix< ValueType > vnlMatrixType;
  typedef vnl_vector< ValueType > vnlVectorType;
  typedef std::vector< std:: vector< std::vector< std::string > > > StringArray3DType;
  typedef Directory DirectoryType;
  typedef typename DirectoryType::Pointer DirectoryPointer;

  typedef TInputImage ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef ImageFileReader< ImageType > ReaderType;
  typedef typename ReaderType::Pointer ReaderPointer;
  typedef itk::ImageRegionIterator< ImageType > IteratorType;
  typedef PermuteAxesImageFilter< ImageType > PermuteAxesFilterType;
  typedef typename PermuteAxesFilterType::Pointer PermuteAxesFilterPointer;


  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::RegionType RegionType;
  typedef typename ImageType::PointType PointType;
  typedef typename SizeType::SizeValueType SizeValueType;

  void Read( std::istream& os );
  void UpdateFileNameLookup( std::istream& os );
  void CreateStitchImage();
  void AllocateROI();

  void GetSettingFieldName( StringVectorType& name )
  {
    name = m_SettingName;
  }

  void GetSettingFieldValue( DoubleVectorType& name )
  {
    name = m_SettingValue;
  }

  unsigned int GetTotalNumberOfTiles()
  {
    return this->m_NumberOfTiles;
  }

  unsigned int * GetTileNumber()
  {
    return m_TileNumber;
  }

  ValueType * GetTileSize()
  {
    return  m_TileSize;
  }

  PointType & GetMinimumStart()
  {
    return m_MinimumStart;
  }

  PointType & GetMaximumEnd()
  {
    return m_MaximumEnd;
  }

  SizeType& GetTileDimension()
  {
    return m_TileDimension;
  }

  SpacingType& GetTileSpacing()
  {
    return m_TileSpacing;
  }

  void SetTileDirectory( std::string iDirectory )
  {
    m_Directory = iDirectory;
  }


  void SetChannelNumber( unsigned int iChannel )
  {
    m_ChannelNumber = iChannel;
  }

  void SetTimePoint( unsigned int iTp )
  {
    m_TimePoint = iTp;
  }

  ValueType * GetStitchSize()
  {
    return  m_StitchSize;
  }

  PointType & GetStitchOrigin()
  {
    return m_StitchOrigin;
  }

  SizeType & GetStitchDimension()
  {
    return m_StitchDimension;
  }

  IndexType & GetStitchIndex()
  {
    return m_StitchIndex;
  }

  RegionType & GetStitchRegion()
  {
    return m_StitchRegion;
  }

  ImagePointer GetStitchImage()
  {
    return m_StitchedImage;
  }

  ImagePointer GetROIImage()
  {
    return m_ROIImage;
  }

  void SetROIOrigin( PointType& iOrigin )
  {
    m_ROIOrigin = iOrigin;
  }

  void SetROI( RegionType& iRegion )
  {
    m_ROI = iRegion;
  }

  StringArray3DType& GetTileFileNameArray()
  {
    return m_TileFileNameArray;
  }

  std::string m_Path;
  std::string m_Directory;
  std::string m_ChannelName;
  unsigned int m_ChannelNumber;
  unsigned int m_TimePoint;

  protected:
  SettingsInfoExtractionFilter();
  ~SettingsInfoExtractionFilter() {}

  void UpdateStitchDimensions( std::istream& os );
  void UpdateTileCoverage( std::istream& os );
  void TransformCoordinateAxes();
  void ReadTileInfo( std::istream& os );
  void FillROI();
  void OverlapRegion( ImagePointer A, ImagePointer B,
                      RegionType& rA, RegionType& rB );

  unsigned int m_Dimension;
  StringVectorType m_SettingName;
  DoubleVectorType m_SettingValue;

  unsigned int m_NumberOfTiles;
  unsigned int m_TileNumber[3];
  ValueType m_TileSize[3];
  SizeType m_TileDimension;
  SpacingType m_TileSpacing;
  StringArray3DType m_TileFileNameArray;

  PointType m_MinimumStart;
  PointType m_MaximumEnd;

  DoubleVectorType m_TileCoverStart[3];
  DoubleVectorType m_TileCoverEnd[3];
  unsigned int m_ScanStart[3];
  unsigned int m_ScanEnd[3];

  vnlMatrixType m_TileInfoValue;
  vnlMatrixType m_TransformedTileInfoValue;

  ImagePointer m_StitchedImage;
  double     m_StitchSize[3];
  PointType  m_StitchOrigin;
  SizeType   m_StitchDimension;
  IndexType  m_StitchIndex;
  RegionType m_StitchRegion;

  ImagePointer m_ROIImage;
  PointType m_ROIOrigin;
  RegionType m_ROI;

  private:
    SettingsInfoExtractionFilter ( Self& );   // intentionally not implemented
    void operator= ( const Self& );   // intentionally not implemented
  };

} /* namespace itk */

#include "itkSettingsInfoExtractionFilter.txx"
#endif
