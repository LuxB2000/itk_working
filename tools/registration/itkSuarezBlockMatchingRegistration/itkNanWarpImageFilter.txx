/*=========================================================================

  This work has been made to be used with the libraries Insight Segmentation
  & Registration Toolkit 2.8.1.

  Program:    itkSuarezBlockMatchingRegistration
  Module:     itkVectorDivideImageFilter.txx

  This work has been

  Derived from: Program: Insight Segmentation & Registration Toolkit
                Module:  $RCSfile: itkNanWarpImageFilter.txx,v $
                Date:    $Date: 2005/07/27 15:21:12 $
                Version: $Revision: 1.25 $
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

#ifndef __itkNanWarpImageFilter_txx
#define __itkNanWarpImageFilter_txx
#include "itkNanWarpImageFilter.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkNumericTraits.h"
#include "itkProgressReporter.h"

namespace itk
{

/**
 * Default constructor.
 */
template <class TInputImage,class TOutputImage,class TDeformationField>
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::NanWarpImageFilter()
{
  // Setup the number of required inputs
  this->SetNumberOfRequiredInputs( 2 );  
  
  // Setup default values
  m_OutputSpacing.Fill( 1.0 );
  m_OutputOrigin.Fill( 0.0 );

  m_EdgePaddingValue = NumericTraits<PixelType>::Zero;

  // Setup default interpolator
  typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();

  m_Interpolator = 
    static_cast<InterpolatorType*>( interp.GetPointer() );

}

/**
 * Standard PrintSelf method.
 */
template <class TInputImage,class TOutputImage,class TDeformationField>
void
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::PrintSelf(std::ostream& os, Indent indent) const
{

  Superclass::PrintSelf(os, indent);

  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "EdgePaddingValue: "
     << static_cast<typename NumericTraits<PixelType>::PrintType>(m_EdgePaddingValue)
     << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  
}


/**
 * Set the output image spacing.
 *
 */
template <class TInputImage,class TOutputImage,class TDeformationField>
void
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::SetOutputSpacing(
  const double* spacing)
{
  SpacingType s(spacing);
  this->SetOutputSpacing( s );
}


/**
 * Set the output image origin.
 *
 */
template <class TInputImage,class TOutputImage,class TDeformationField>
void
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::SetOutputOrigin(
  const double* origin)
{
  PointType p(origin);
  this->SetOutputOrigin(p);
}



/**
 * Set deformation field as Inputs[1] for this ProcessObject.
 *
 */
template <class TInputImage,class TOutputImage,class TDeformationField>
void
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::SetDeformationField(
  const DeformationFieldType * field )
{
  // const cast is needed because the pipeline is not const-correct.
  DeformationFieldType * input =  
       const_cast< DeformationFieldType * >( field );
  this->ProcessObject::SetNthInput( 1, input );
}


/**
 * Return a pointer to the deformation field.
 */
template <class TInputImage,class TOutputImage,class TDeformationField>
typename NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::DeformationFieldType *
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::GetDeformationField(void)
{
  return static_cast<DeformationFieldType *>
    ( this->ProcessObject::GetInput( 1 ));
}


/**
 * Setup state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be setup before ThreadedGenerateData
 */
template <class TInputImage,class TOutputImage,class TDeformationField>
void
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::BeforeThreadedGenerateData()
{

  if( !m_Interpolator )
    {
    itkExceptionMacro(<< "Interpolator not set");
    }

  // Connect input image to interpolator
  m_Interpolator->SetInputImage( this->GetInput() );

}

/**
 * Setup state of filter after multi-threading.
 */
template <class TInputImage,class TOutputImage,class TDeformationField>
void
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::AfterThreadedGenerateData()
{
  // Disconnect input image from interpolator
  m_Interpolator->SetInputImage( NULL );

}


/**
 * Compute the output for the region specified by outputRegionForThread.
 */
template <class TInputImage,class TOutputImage,class TDeformationField>
void
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::ThreadedGenerateData(
  const OutputImageRegionType& outputRegionForThread,
  int threadId )
{

  InputImageConstPointer inputPtr = this->GetInput();
  OutputImagePointer outputPtr = this->GetOutput();
  DeformationFieldPointer fieldPtr = this->GetDeformationField();

  // support progress methods/callbacks
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  
  // iterator for the output image
  ImageRegionIteratorWithIndex<OutputImageType> outputIt(
    outputPtr, outputRegionForThread );

  // iterator for the deformation field
  ImageRegionIterator<DeformationFieldType> fieldIt(
    fieldPtr, outputRegionForThread );
  
  IndexType index;
  PointType point;
  DisplacementType displacement;

  while( !outputIt.IsAtEnd() )
    {
    // get the required displacement
    displacement = fieldIt.Get();

    // if is real displacement (not nan)
    if ( displacement[0] == displacement[0] )
      {

      // get the output image index
      index = outputIt.GetIndex();
      outputPtr->TransformIndexToPhysicalPoint( index, point );
      
      // compute the required input image point
      for(unsigned int j = 0; j < ImageDimension; j++ )
	{
	point[j] += displacement[j];
	}
      
      // get the interpolated value
      if( m_Interpolator->IsInsideBuffer( point ) )
	{
	PixelType value = static_cast<PixelType>( 
	  m_Interpolator->Evaluate( point ) );
	outputIt.Set( value );
	}
      else
	{
	outputIt.Set( m_EdgePaddingValue );
	}
      }
    else
      {
      outputIt.Set( m_EdgePaddingValue );
      }
      ++outputIt;
      ++fieldIt; 
      progress.CompletedPixel();
      }

}


template <class TInputImage,class TOutputImage,class TDeformationField>
void
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::GenerateInputRequestedRegion()
{

  // call the superclass's implementation
  Superclass::GenerateInputRequestedRegion();

  // request the largest possible region for the input image
  InputImagePointer inputPtr = 
    const_cast< InputImageType * >( this->GetInput() );

  if( inputPtr )
    {
    inputPtr->SetRequestedRegionToLargestPossibleRegion();
    }

  // just propagate up the output requested region for the 
  // deformation field.
  DeformationFieldPointer fieldPtr = this->GetDeformationField();
  OutputImagePointer outputPtr = this->GetOutput();
  if( fieldPtr )
    {
    fieldPtr->SetRequestedRegion( outputPtr->GetRequestedRegion() );
    }

}


template <class TInputImage,class TOutputImage,class TDeformationField>
void
NanWarpImageFilter<TInputImage,TOutputImage,TDeformationField>
::GenerateOutputInformation()
{
  // call the superclass's implementation of this method
  Superclass::GenerateOutputInformation();

  OutputImagePointer outputPtr = this->GetOutput();

  outputPtr->SetSpacing( m_OutputSpacing );
  outputPtr->SetOrigin( m_OutputOrigin );

  DeformationFieldPointer fieldPtr = this->GetDeformationField();
  if( fieldPtr )
    {
    outputPtr->SetLargestPossibleRegion( fieldPtr->
                                         GetLargestPossibleRegion() );
    }

}


} // end namespace itk

#endif
 
