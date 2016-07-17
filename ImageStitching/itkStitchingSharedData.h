#ifndef __itkStitchingSharedData_h_
#define __itkStitchingSharedData_h_

#include <vector>
#include <cstring>
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
  typedef ImageFileReader< ImageType > ReaderType;
  typedef typename ReaderType::Pointer ReaderPointer;

  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::SpacingType SpacingType;

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

  itkSetMacro( DeconvolutionIterations,  unsigned int );
  itkSetMacro( TileDimension,            SizeType );
  itkSetMacro( TileSpacing,              SpacingType );
  itkSetMacro( PSFPath,                  std::string );
  itkGetConstMacro( TileFileNameArray,   StringArray3DType );
  itkGetConstMacro( CorrectionThreshold, double );

  void SetCorrectionInfo( std::string& iName, double& var, double& thresh )
  {
    m_CorrectionFilename = iName;
    m_CorrectionVariance = var;
    m_CorrectionThreshold = thresh;

    ReadCorrectionImage();
  }

  void ReadPSFImage();

  StringArray3DType m_TileFileNameArray;

  RImagePointer m_CorrectionImage;
  double        m_CorrectionThreshold;
  std::string   m_CorrectionFilename;
  double        m_CorrectionVariance;

  ImagePointer  m_PSF;
  std::string   m_PSFPath;
  unsigned int  m_DeconvolutionIterations;

  DoubleVectorType  m_TileCover[3][2][2];
  // Dimension, Start/End, Clipping

protected:
  StitchingSharedData();
  ~StitchingSharedData(){}
  void ReadCorrectionImage();
  void PrintSelf(std::ostream& os, Indent indent) const;

  SizeType          m_TileDimension;
  SpacingType       m_TileSpacing;

private:
  StitchingSharedData(const Self&) {}
  void operator=(const Self&) {}
};

} // end namespace itk

#include "itkStitchingSharedData.txx"
#endif

