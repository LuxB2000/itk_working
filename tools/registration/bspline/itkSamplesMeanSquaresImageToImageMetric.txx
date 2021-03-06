/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSamplesMeanSquaresImageToImageMetric.txx,v $
  Language:  C++
  Date:      $Date: 2004/12/21 22:47:27 $
  Version:   $Revision: 1.46 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkSamplesMeanSquaresImageToImageMetric_txx
#define _itkSamplesMeanSquaresImageToImageMetric_txx

#include "itkSamplesMeanSquaresImageToImageMetric.h"
#include "itkImageRandomConstIteratorWithIndex.h"

namespace itk
{

unsigned long m_MS_NumberOfSpatialSamples;

/*
 * Constructor
 */
template <class TFixedImage, class TMovingImage> 
SamplesMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::SamplesMeanSquaresImageToImageMetric()
{
  itkDebugMacro("Constructor");
  m_MS_NumberOfSpatialSamples = 1000;
}

/*
 * Get the match Measure
 */
template <class TFixedImage, class TMovingImage> 
typename SamplesMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>::MeasureType
SamplesMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValue( const TransformParametersType & parameters ) const
{

  itkDebugMacro("GetValue( " << parameters << " ) ");

  FixedImageConstPointer fixedImage = this->m_FixedImage;

  if( !fixedImage ) 
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  typedef  itk::ImageRandomConstIteratorWithIndex<FixedImageType> FixedIteratorType;

  FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );
  ti.SetNumberOfSamples(m_MS_NumberOfSpatialSamples);

  typename FixedImageType::IndexType index;

  MeasureType measure = NumericTraits< MeasureType >::Zero;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  while(!ti.IsAtEnd())
    {

    index = ti.GetIndex();
    
    typename Superclass::InputPointType inputPoint;
    fixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

    if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
      {
      ++ti;
      continue;
      }

    typename Superclass::OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside( transformedPoint ) )
      {
      ++ti;
      continue;
      }

    if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
      {
      const RealType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );
      const RealType fixedValue   = ti.Get();
      this->m_NumberOfPixelsCounted++;
      const RealType diff = movingValue - fixedValue; 
      measure += diff * diff; 
      }

    ++ti;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<<"All the points mapped to outside of the moving image");
    }
  else
    {
    measure /= this->m_NumberOfPixelsCounted;
    }

  //std::cout << measure << std::endl;

  return measure;

}





/*
 * Get the Derivative Measure
 */
template < class TFixedImage, class TMovingImage> 
void
SamplesMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative  ) const
{

  itkDebugMacro("GetDerivative( " << parameters << " ) ");
  
  if( !this->GetGradientImage() )
    {
    itkExceptionMacro(<<"The gradient image is null, maybe you forgot to call Initialize()");
    }

  FixedImageConstPointer fixedImage = this->m_FixedImage;

  if( !fixedImage ) 
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  const unsigned int ImageDimension = FixedImageType::ImageDimension;


  typedef  itk::ImageRandomConstIteratorWithIndex<FixedImageType> FixedIteratorType;

  FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );
  ti.SetNumberOfSamples(m_MS_NumberOfSpatialSamples);

//   typedef  itk::ImageRegionConstIteratorWithIndex<
//     ITK_TYPENAME Superclass::GradientImageType> GradientIteratorType;

  typename FixedImageType::IndexType index;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType( ParametersDimension );
  derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  ti.GoToBegin();

  while(!ti.IsAtEnd())
    {

    index = ti.GetIndex();
    
    typename Superclass::InputPointType inputPoint;
    fixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

    if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
      {
      ++ti;
      continue;
      }

    typename Superclass::OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside( transformedPoint ) )
      {
      ++ti;
      continue;
      }

    if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
      {
      const RealType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );

      const TransformJacobianType & jacobian =
        this->m_Transform->GetJacobian( inputPoint ); 

      
      const RealType fixedValue     = ti.Value();
      this->m_NumberOfPixelsCounted++;
      const RealType diff = movingValue - fixedValue; 

      // Get the gradient by NearestNeighboorInterpolation: 
      // which is equivalent to round up the point components.
      typedef typename Superclass::OutputPointType OutputPointType;
      typedef typename OutputPointType::CoordRepType CoordRepType;
      typedef ContinuousIndex<CoordRepType,MovingImageType::ImageDimension>
        MovingImageContinuousIndexType;

      MovingImageContinuousIndexType tempIndex;
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex( transformedPoint, tempIndex );

      typename MovingImageType::IndexType mappedIndex; 
      for( unsigned int j = 0; j < MovingImageType::ImageDimension; j++ )
        {
        mappedIndex[j] = static_cast<long>( vnl_math_rnd( tempIndex[j] ) );
        }

      const GradientPixelType gradient = 
        this->GetGradientImage()->GetPixel( mappedIndex );

      for(unsigned int par=0; par<ParametersDimension; par++)
        {
        RealType sum = NumericTraits< RealType >::Zero;
        for(unsigned int dim=0; dim<ImageDimension; dim++)
          {
          sum += 2.0 * diff * jacobian( dim, par ) * gradient[dim];
          }
        derivative[par] += sum;
        }
      }

    ++ti;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<<"All the points mapped to outside of the moving image");
    }
  else
    {
    for(unsigned int i=0; i<ParametersDimension; i++)
      {
      derivative[i] /= this->m_NumberOfPixelsCounted;
      }
    }

}



// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  > void
SamplesMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>::ReinitializeSeed()
{
  // This method should be the same used in the ImageRandomIterator
  //vnl_sample_reseed();
Statistics::MersenneTwisterRandomVariateGenerator::GetInstance()->SetSeed();
}

// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  > void
SamplesMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::ReinitializeSeed(int seed)
{
  // This method should be the same used in the ImageRandomIterator
  //vnl_sample_reseed(seed);
Statistics::MersenneTwisterRandomVariateGenerator::GetInstance()->SetSeed(seed);
}




/*
 * Set the number of samples
 */
template <class TFixedImage, class TMovingImage> 
void
SamplesMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::SetNumberOfSamples(unsigned long num) const 
{
m_MS_NumberOfSpatialSamples = num;
}





/*
 * Get both the match Measure and theDerivative Measure 
 */
template <class TFixedImage, class TMovingImage> 
void
SamplesMeanSquaresImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivative(const TransformParametersType & parameters, 
                        MeasureType & value, DerivativeType  & derivative) const
{

  itkDebugMacro("GetValueAndDerivative( " << parameters << " ) ");

  if( !this->GetGradientImage() )
    {
    itkExceptionMacro(<<"The gradient image is null, maybe you forgot to call Initialize()");
    }

  FixedImageConstPointer fixedImage = this->m_FixedImage;

  if( !fixedImage ) 
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  const unsigned int ImageDimension = FixedImageType::ImageDimension;

  typedef  itk::ImageRandomConstIteratorWithIndex<FixedImageType> FixedIteratorType;

  FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );
  
  ti.SetNumberOfSamples(m_MS_NumberOfSpatialSamples);

  typename FixedImageType::IndexType index;

  MeasureType measure = NumericTraits< MeasureType >::Zero;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType( ParametersDimension );
  derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  ti.GoToBegin();

  while(!ti.IsAtEnd())
    {

    index = ti.GetIndex();
    
    typename Superclass::InputPointType inputPoint;
    fixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

    if( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
      {
      ++ti;
      continue;
      }

    typename Superclass::OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    if( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside( transformedPoint ) )
      {
      ++ti;
      continue;
      }

    if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
      {
      const RealType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );

      const TransformJacobianType & jacobian =
        this->m_Transform->GetJacobian( inputPoint ); 

      
      const RealType fixedValue     = ti.Value();
      this->m_NumberOfPixelsCounted++;

      const RealType diff = movingValue - fixedValue; 
  
      measure += diff * diff;

      // Get the gradient by NearestNeighboorInterpolation: 
      // which is equivalent to round up the point components.
      typedef typename Superclass::OutputPointType OutputPointType;
      typedef typename OutputPointType::CoordRepType CoordRepType;
      typedef ContinuousIndex<CoordRepType,MovingImageType::ImageDimension>
        MovingImageContinuousIndexType;

      MovingImageContinuousIndexType tempIndex;
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex( transformedPoint, tempIndex );

      typename MovingImageType::IndexType mappedIndex; 
      for( unsigned int j = 0; j < MovingImageType::ImageDimension; j++ )
        {
        mappedIndex[j] = static_cast<long>( vnl_math_rnd( tempIndex[j] ) );
        }

      const GradientPixelType gradient = 
        this->GetGradientImage()->GetPixel( mappedIndex );

      for(unsigned int par=0; par<ParametersDimension; par++)
        {
        RealType sum = NumericTraits< RealType >::Zero;
        for(unsigned int dim=0; dim<ImageDimension; dim++)
          {
          sum += 2.0 * diff * jacobian( dim, par ) * gradient[dim];
          }
        derivative[par] += sum;
        }
      }

    ++ti;
    }

  if( !this->m_NumberOfPixelsCounted )
    {
    itkExceptionMacro(<<"All the points mapped to outside of the moving image");
    }
  else
    {
    for(unsigned int i=0; i<ParametersDimension; i++)
      {
      derivative[i] /= this->m_NumberOfPixelsCounted;
      }
    measure /= this->m_NumberOfPixelsCounted;
    }

  value = measure;
  
  std::cout << "Mean square = " << value << std::endl;  

}

} // end namespace itk


#endif
