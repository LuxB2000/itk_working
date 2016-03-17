/*=========================================================================

  This work has been made to be used with the libraries Insight Segmentation
  & Registration Toolkit 2.8.1.

  Program:    itkSuarezBlockMatchingRegistration
  Module:     itkNanWarpImageFilter.h

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

#ifndef __itkNanWarpImageFilter_h
#define __itkNanWarpImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkPoint.h"
#include "itkFixedArray.h"

namespace itk
{

/** \class NanWarpImageFilter
 * \brief Warps an image using an input deformation field.
 *
 * NanWarpImageFilter warps an existing image with respect to
 * a given deformation field.
 *
 * A deformation field is represented as a image whose pixel type is some
 * vector type with at least N elements, where N is the dimension of
 * the input image. The vector type must support element access via operator
 * [].
 *
 * The output image is produced by inverse mapping: the output pixels
 * are mapped back onto the input image. This scheme avoids the creation of
 * any holes and overlaps in the output image.
 *
 * Each vector in the deformation field represent the distance between
 * a geometric point in the input space and a point in the output space such 
 * that:
 *
 * \f[ p_{in} = p_{out} + d \f]
 *
 * Typically the mapped position does not correspond to an integer pixel 
 * position in the input image. Interpolation via an image function
 * is used to compute values at non-integer positions. The default 
 * interpolation typed used is the LinearInterpolateImageFunction.
 * The user can specify a particular interpolation function via
 * SetInterpolator(). Note that the input interpolator must derive
 * from base class InterpolateImageFunction.
 *
 * Position mapped to outside of the input image buffer are assigned
 * a edge padding value.
 *
 * The LargetPossibleRegion for the output is inherited from the 
 * input deformation field. The output image spacing and origin may be set 
 * via SetOutputSpacing, SetOutputOrigin. The default are respectively a 
 * vector of 1's and a vector of 0's.
 *
 * This class is templated over the type of the input image, the
 * type of the output image and the type of the deformation field.
 *
 * The input image is set via SetInput. The input deformation field
 * is set via SetDeformationField.
 *
 * This filter is implemented as a multithreaded filter.
 *
 * \warning This filter assumes that the input type, output type
 * and deformation field type all have the same number of dimensions.
 *
 * \ingroup GeometricTransforms MultiThreaded
 */

/* Difference with itkWarpImageFilter
 *
 * There is a modification in the way the warping in done, in the
 * sense that for nan vectors, the warping assumes their edge padding
 * value.
 */

template <
  class TInputImage,
  class TOutputImage,
  class TDeformationField
  >
class ITK_EXPORT NanWarpImageFilter :
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef NanWarpImageFilter      Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkTypeMacro( NanWarpImageFilter, ImageToImageFilter );

  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Inherit some types from the superclass. */
  typedef typename Superclass::InputImageType        InputImageType;
  typedef typename Superclass::InputImagePointer     InputImagePointer;
  typedef typename Superclass::OutputImageType       OutputImageType;
  typedef typename Superclass::OutputImagePointer    OutputImagePointer;
  typedef typename Superclass::InputImageConstPointer InputImageConstPointer;
  typedef typename OutputImageType::IndexType        IndexType;
  typedef typename OutputImageType::SizeType         SizeType;
  typedef typename OutputImageType::PixelType        PixelType;
  typedef typename OutputImageType::SpacingType      SpacingType;

  /** Determine the image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension );
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension );
  itkStaticConstMacro(DeformationFieldDimension, unsigned int,
                      TDeformationField::ImageDimension );

  /** Deformation field typedef support. */
  typedef TDeformationField    DeformationFieldType;
  typedef typename DeformationFieldType::Pointer  DeformationFieldPointer;
  typedef typename DeformationFieldType::PixelType DisplacementType;

  /** Interpolator typedef support. */
  typedef double CoordRepType;
  typedef InterpolateImageFunction<InputImageType,CoordRepType>   InterpolatorType;
  typedef typename InterpolatorType::Pointer   InterpolatorPointer;
  typedef LinearInterpolateImageFunction<InputImageType,CoordRepType>
  DefaultInterpolatorType;

  /** Point type */
  typedef Point<CoordRepType,itkGetStaticConstMacro(ImageDimension)> PointType;

  /** Set the deformation field. */
  void SetDeformationField( const DeformationFieldType * field );

  /** Get a pointer the deformation field. */
  DeformationFieldType * GetDeformationField(void);

  /** Set the interpolator function. */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the interpolator function. */
  itkGetObjectMacro( Interpolator, InterpolatorType );

  /** Set the output image spacing. */
  itkSetMacro(OutputSpacing, SpacingType);
  virtual void SetOutputSpacing( const double* values);

  /** Get the output image spacing. */
  itkGetConstReferenceMacro(OutputSpacing, SpacingType);

  /** Set the output image origin. */
  itkSetMacro(OutputOrigin, PointType);
  virtual void SetOutputOrigin( const double* values);

  /** Get the output image origin. */
  itkGetConstReferenceMacro(OutputOrigin, PointType);

  /** Set the edge padding value */
  itkSetMacro( EdgePaddingValue, PixelType );

  /** Get the edge padding value */
  itkGetMacro( EdgePaddingValue, PixelType );

  /** NanWarpImageFilter produces an image which is a different
   * size than its input image. As such, it needs to provide an
   * implemenation for GenerateOutputInformation() which set
   * the output information according the OutputSpacing, OutputOrigin
   * and the deformation field's LargestPossibleRegion. */
  virtual void GenerateOutputInformation();

  /** It is difficult to compute in advance the input image region
   * required to compute the requested output region. Thus the safest
   * thing to do is to request for the whole input image.
   *
   * For the deformation field, the input requested region
   * set to be the same as that of the output requested region. */
  virtual void GenerateInputRequestedRegion();

  /** This method is used to set the state of the filter before 
   * multi-threading. */
  virtual void BeforeThreadedGenerateData();

  /** This method is used to set the state of the filter after 
   * multi-threading. */
  virtual void AfterThreadedGenerateData();

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(SameDimensionCheck1,
    (Concept::SameDimension<ImageDimension, InputImageDimension>));
  itkConceptMacro(SameDimensionCheck2,
    (Concept::SameDimension<ImageDimension, DeformationFieldDimension>));
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<typename TInputImage::PixelType>));
  itkConceptMacro(DeformationFieldHasNumericTraitsCheck,
    (Concept::HasNumericTraits<typename TDeformationField::PixelType::ValueType>));
  /** End concept checking */
#endif

protected:
  NanWarpImageFilter();
  ~NanWarpImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** NanWarpImageFilter is implemented as a multi-threaded filter.
   * As such, it needs to provide and implementation for 
   * ThreadedGenerateData(). */
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            int threadId );

private:
  NanWarpImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  PixelType                  m_EdgePaddingValue;
  SpacingType                m_OutputSpacing;
  PointType                  m_OutputOrigin;

  InterpolatorPointer        m_Interpolator;
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNanWarpImageFilter.txx"
#endif

#endif
