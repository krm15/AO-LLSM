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


namespace itk
{
template < class TValueType >
class ITK_EXPORT SettingsInfoExtractionFilter : public LightObject
{
  public:
    typedef SettingsInfoExtractionFilter Self;
    typedef LightObject                  Superclass;
    typedef SmartPointer< Self >         Pointer;
    typedef SmartPointer< const Self >   ConstPointer;

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
    typedef itk::Directory DirectoryType;

    void Read( std::istream& os );

  protected:
  SettingsInfoExtractionFilter();
  ~SettingsInfoExtractionFilter() {}

  unsigned int m_Dimension;
  StringVectorType m_SettingName;
  DoubleVectorType m_SettingValue;

  private:
    SettingsInfoExtractionFilter ( Self& );   // intentionally not implemented
    void operator= ( const Self& );   // intentionally not implemented
  };

} /* namespace itk */

#include "itkSettingsInfoExtractionFilter.txx"
#endif
