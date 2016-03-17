/*=========================================================================

  This work has been made to be used with the libraries Insight Segmentation
  & Registration Toolkit 2.8.1.

  Program:    itkSuarezBlockMatchingRegistration
  Module:     itkVectorDivideImageFilter.h

  This work has been

  Derived from: Program: Insight Segmentation & Registration Toolkit
                Module:  $RCSfile: itkWarpImageFilter.h,v $
                Date:    $Date: 2006/04/04 13:13:51 $
                Version: $Revision: 1.21 $
                License of the original work:

            Copyright (c) Insight Software Consortium. All rights reserved.
            See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
            
            This software is distributed WITHOUT ANY WARRANTY; without even 
            the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
            PURPOSE.  See the above copyright notices for more information.

  Modified by: Eduardo Suarez
               Instituto Tecnologico de Canarias - Gobierno de Canarias
               http://www.itccanarias.org/
  Date:        2006/08/07


  License of the contributions of this derivative work:

    Copyright (c) 2006 Instituto Tecnologico de Canarias - Gobierno de Canarias
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

     * Redistributions of source code must retain the above copyright notice,
       this list of conditions and the following disclaimer.

     * Redistributions in binary form must reproduce the above copyright notice,
       this list of conditions and the following disclaimer in the documentation
       and/or other materials provided with the distribution.

     * The name of the Instituto Tecnologico de Canarias - Gobierno de Canarias,
       may be used to endorse or promote products derived from this software
       without specific prior written permission.
    
     * Modified source versions must be plainly marked as such, and must not be
       misrepresented as being the original software.
    
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


  References:

  [1] Nonrigid Registration using Regularized Matching Weighted by
      Local Structure, E. Suarez, C.-F. Westin, E. Rovaris, J. Ruiz-Alzola
      MICCAI, Tokyo, Japan. 2002. Pages 581-589. 2002

  [2] Fast Entropy-Based Nonrigid Registration. E. Suarez, J. A. Santana,
      E. Rovaris, C.-F. Westin, J. Ruiz-Alzola. Computer Aided Systems Theory 
      (EUROCAST'03), Lecture Notes in Computer Science 2809.
      Las Palmas de Gran Canaria, Spain. Pages 607-615. 2003


  Acknowledgements:

    This work has been partially funded by project Torres Quevedo
    PTQ2004-1444 of the spanish goverment.

    The authors thanks to the National Alliance of Medical Image
    Computing (http://www.na-mic.org/) for his technical support in 
    the Insight Segmentation & Registration Toolkit.

    The authors also thank to Dan Blezek for his first approach to an
    implementation of this filter, to Luis Ibanez for his technical
    support.

=========================================================================*/

#ifndef __itkVectorDivideImageFilter_h
#define __itkVectorDivideImageFilter_h

#include "itkBinaryFunctorImageFilter.h"
#include "itkNumericTraits.h"

namespace itk
{
  
/** \class DivideImageFilter
 * \brief Implements an operator for pixel-wise division of two images.
 *
 * This class is parametrized over the types of the two 
 * input images and the type of the output image. When the divisor is zero,
 * the division result is set to the maximum number that can be represneted  by default to 
 * avoid exception. Numeric conversions (castings) are done by the C++ defaults.
 * 
 * \ingroup IntensityImageFilters  Multithreaded
 */

namespace Function {  
  
template< class TInput1, class TInput2, class TOutput>
class VectorDiv
{
public:
  VectorDiv() {};
  ~VectorDiv() {};
  bool operator!=( const VectorDiv & ) const
  {
    return false;
  }
  bool operator==( const VectorDiv & other ) const
  {
    return !(*this != other);
  }
  inline TOutput operator()( const TInput1 & A, const TInput2 & B)
  {
    return (TOutput)(A / B);
//     if(B != (TInput2) 0)
//       return (TOutput)(A / B);
//     else
//       return NumericTraits<TOutput>::max();
  }
}; 
}

template <class TInputImage1, class TInputImage2, class TOutputImage>
class ITK_EXPORT VectorDivideImageFilter :
    public
BinaryFunctorImageFilter<TInputImage1,TInputImage2,TOutputImage, 
                         Function::VectorDiv< 
  typename TInputImage1::PixelType, 
  typename TInputImage2::PixelType,
  typename TOutputImage::PixelType>   >
{
public:
  /**
   * Standard "Self" typedef.
   */
  typedef VectorDivideImageFilter  Self;

  /**
   * Standard "Superclass" typedef.
   */
  typedef BinaryFunctorImageFilter<TInputImage1,TInputImage2,TOutputImage, 
                                   Function::VectorDiv< 
    typename TInputImage1::PixelType, 
    typename TInputImage2::PixelType,
    typename TOutputImage::PixelType>   
  > Superclass;

  /** 
   * Smart pointer typedef support 
   */
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /**
   * Method for creation through the object factory.
   */
  itkNewMacro(Self);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(IntConvertibleToInput2Check,
    (Concept::Convertible<int, typename TInputImage2::PixelType>));
  itkConceptMacro(Input1Input2OutputVectorDivisionOperatorsCheck,
    (Concept::VectorDivisionOperators<typename TInputImage1::PixelType,
                                typename TInputImage2::PixelType,
                                typename TOutputImage::PixelType>));
  /** End concept checking */
#endif

protected:
  VectorDivideImageFilter() {}
  virtual ~VectorDivideImageFilter() {}
  VectorDivideImageFilter(const Self&) {}
  void operator=(const Self&) {}

};

} // end namespace itk


#endif
