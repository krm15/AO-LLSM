#ifndef __itkStitchingSharedData_h_
#define __itkStitchingSharedData_h_

#include <vector>
#include <cstring>
#include "itkDirectory.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

namespace itk
{
	
template <class TInputImage>
class StitchingSharedData : public Object
{
 public:
  typedef StitchingSharedData               Self;
  typedef Object                            Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  itkStaticConstMacro ( ImageDimension, unsigned int,
    TInputImage::ImageDimension );

  /** Run-time type information (and related methods). */
  itkTypeMacro( StitchingSharedData, Object );

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  typedef std::vector< std:: vector< std::vector< std::string > > > StringArray3DType;
  typedef std::vector< double > DoubleVectorType;
  typedef TInputImage ImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::PixelType PixelType;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::RegionType RegionType;
  typedef typename ImageType::PointType PointType;
  typedef typename SizeType::SizeValueType SizeValueType;

  typedef ImageFileReader< ImageType > ReaderType;
  typedef typename ReaderType::Pointer ReaderPointer;
  typedef ImageRegionIterator< ImageType > IteratorType;

  typedef Image< double, 2 > RImageType;
  typedef typename RImageType::Pointer RImagePointer;
  typedef ImageFileReader< RImageType > RReaderType;
  typedef typename RReaderType::Pointer RReaderPointer;
  typedef ImageRegionIterator< RImageType > RIteratorType;

  typedef typename RImageType::SizeType RSizeType;
  typedef typename RImageType::SpacingType RSpacingType;
  typedef typename RImageType::IndexType RIndexType;
  typedef typename RImageType::RegionType RRegionType;
  typedef typename RImageType::PointType RPointType;

  typedef DiscreteGaussianImageFilter< RImageType, RImageType > GaussianFilterType;
  typedef typename GaussianFilterType::Pointer GaussianFilterPointer;
  typedef RegionOfInterestImageFilter< RImageType, RImageType > ROIFilterType;
  typedef typename ROIFilterType::Pointer ROIFilterPointer;
  typedef RescaleIntensityImageFilter< RImageType, RImageType > RescaleFilterType;
  typedef typename RescaleFilterType::Pointer RescaleFilterPointer;

  typedef Directory DirectoryType;
  typedef typename DirectoryType::Pointer DirectoryPointer;

  itkSetMacro( DeconvolutionIterations,  unsigned int );
  itkSetMacro( TileDimension,            SizeType );
  itkSetMacro( TileSpacing,              SpacingType );
  itkSetMacro( PSFPath,                  std::string );
  itkGetConstMacro( TileFileNameArray,   StringArray3DType );
  itkGetConstMacro( CorrectionThreshold, double );
  itkGetConstMacro( CorrectionVariance, double );
  itkGetConstMacro( CorrectionFilename, std::string );

  itkGetConstMacro( TileDimension,      SizeType );
  itkGetConstMacro( TileOverlap,        PointType );
  itkGetConstMacro( TileSize,           PointType );
  itkGetConstMacro( TileNumber,         IndexType );
  itkGetConstMacro( TileSpacing,        SpacingType );


  void SetCorrectionInfo( std::string& iName, double& var, double& thresh )
  {
    m_CorrectionFilename = iName;
    m_CorrectionVariance = var;
    m_CorrectionThreshold = thresh;

    ReadCorrectionImage();
  }

  void ReadPSFImage();
  void ComputeScalingFactor( std::string& iDirName, std::string& iSearchStringCH );

  StringArray3DType m_TileFileNameArray;

  RImagePointer m_CorrectionImage;
  double        m_CorrectionThreshold;
  std::string   m_CorrectionFilename;
  double        m_CorrectionVariance;

  double        m_ScalingFactor;

  IndexType         m_TileNumber;
  PointType         m_TileSize;
  SizeType          m_TileDimension;
  SpacingType       m_TileSpacing;
  PointType         m_TileOverlap;

  ImagePointer  m_PSF;
  std::string   m_PSFPath;
  unsigned int  m_DeconvolutionIterations;

  DoubleVectorType  m_TileCover[3][2][2];
  // Dimension, Start/End, Clipping
  DoubleVectorType  m_TileOffset[3];
  DoubleVectorType  m_TileEffectiveOffset[3];

protected:
  StitchingSharedData();
  ~StitchingSharedData(){}
  void ReadCorrectionImage();
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  StitchingSharedData(const Self&) {}
  void operator=(const Self&) {}
};

} // end namespace itk

#include "itkStitchingSharedData.txx"
#endif

