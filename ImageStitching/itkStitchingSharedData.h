#ifndef __itkStitchingSharedData_h_
#define __itkStitchingSharedData_h_

#include <vector>
#include <cstring>
#include "itkImage.h"

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
  typedef Image< double, 2 > RImageType;
  typedef typename RImageType::Pointer RImagePointer;

  itkGetConstMacro( TileFileNameArray,  StringArray3DType );

  StringArray3DType m_TileFileNameArray;
  DoubleVectorType  m_TileCover[3][2][2];
  RImagePointer m_CorrectionImage;
  double        m_CorrectionThreshold;
  ImagePointer  m_PSF;

  // Dimension, Start/End, Clipping

protected:
  StitchingSharedData()
  {
    m_CorrectionThreshold = 1200;
    m_CorrectionImage = ITK_NULLPTR;
    m_PSF = ITK_NULLPTR;
  }


  ~StitchingSharedData()
  {
	}


private:
  StitchingSharedData(const Self&) {}
  void operator=(const Self&) {}
};

} // end namespace itk

#endif

