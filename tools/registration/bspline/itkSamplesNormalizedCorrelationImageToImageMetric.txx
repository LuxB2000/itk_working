/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkSamplesNormalizedCorrelationImageToImageMetric.txx,v $
  Language:  C++
  Date:      $Date: 2005/08/10 00:11:24 $
  Version:   $Revision: 1.45 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkSamplesNormalizedCorrelationImageToImageMetric_txx
#define _itkSamplesNormalizedCorrelationImageToImageMetric_txx

#include "itkSamplesNormalizedCorrelationImageToImageMetric.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRandomConstIteratorWithIndex.h"

namespace itk
{

unsigned long m_NumberOfSpatialSamples;

/*
 * Constructor
 */
template <class TFixedImage, class TMovingImage> 
SamplesNormalizedCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::SamplesNormalizedCorrelationImageToImageMetric()
{
  m_SubtractMean = false;
  m_NumberOfSpatialSamples = 100000;
}

/*
 * Get the match Measure
 */
template <class TFixedImage, class TMovingImage> 
typename SamplesNormalizedCorrelationImageToImageMetric<TFixedImage,TMovingImage>::MeasureType
SamplesNormalizedCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::GetValue( const TransformParametersType & parameters ) const
{

  FixedImageConstPointer fixedImage = this->m_FixedImage;

  if( !fixedImage ) 
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  typedef  itk::ImageRandomConstIteratorWithIndex<FixedImageType> FixedIteratorType;


  FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );
  ti.SetNumberOfSamples(m_NumberOfSpatialSamples);

  typename FixedImageType::IndexType index;

  MeasureType measure;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  typedef  typename NumericTraits< MeasureType >::AccumulateType AccumulateType;

  AccumulateType sff = NumericTraits< AccumulateType >::Zero;
  AccumulateType smm = NumericTraits< AccumulateType >::Zero;
  AccumulateType sfm = NumericTraits< AccumulateType >::Zero;
  AccumulateType sf  = NumericTraits< AccumulateType >::Zero;
  AccumulateType sm  = NumericTraits< AccumulateType >::Zero;

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
      sff += fixedValue  * fixedValue;
      smm += movingValue * movingValue;
      sfm += fixedValue  * movingValue;
      if ( this->m_SubtractMean )
        {
        sf += fixedValue;
        sm += movingValue;
        }
      this->m_NumberOfPixelsCounted++;
      }

    ++ti;
    }

  if ( this->m_SubtractMean && this->m_NumberOfPixelsCounted > 0 )
    {
    sff -= ( sf * sf / this->m_NumberOfPixelsCounted );
    smm -= ( sm * sm / this->m_NumberOfPixelsCounted );
    sfm -= ( sf * sm / this->m_NumberOfPixelsCounted );
    }

  const RealType denom = -1.0 * vcl_sqrt(sff * smm );

  if( this->m_NumberOfPixelsCounted > 0 && denom != 0.0)
    {
    measure = sfm / denom;
    }
  else
    {
    measure = NumericTraits< MeasureType >::Zero;
    }

  std::cout << "Normalized correlation = " << -measure << std::endl;  
  
  return measure;

}





/*
 * Get the Derivative Measure
 */
template < class TFixedImage, class TMovingImage> 
void
SamplesNormalizedCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::GetDerivative( const TransformParametersType & parameters,
                 DerivativeType & derivative ) const
{

  if( !this->GetGradientImage() )
    {
    itkExceptionMacro(<<"The gradient image is null, maybe you forgot to call Initialize()");
    }

  FixedImageConstPointer fixedImage = this->m_FixedImage;

  if( !fixedImage ) 
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  const unsigned int dimension = FixedImageType::ImageDimension;

  typedef  itk::ImageRandomConstIteratorWithIndex<FixedImageType> FixedIteratorType;

  FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );
  ti.SetNumberOfSamples(m_NumberOfSpatialSamples);

  typename FixedImageType::IndexType index;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  typedef  typename NumericTraits< MeasureType >::AccumulateType AccumulateType;

  AccumulateType sff  = NumericTraits< AccumulateType >::Zero;
  AccumulateType smm  = NumericTraits< AccumulateType >::Zero;
  AccumulateType sfm = NumericTraits< AccumulateType >::Zero;
  AccumulateType sf  = NumericTraits< AccumulateType >::Zero;
  AccumulateType sm  = NumericTraits< AccumulateType >::Zero;

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType( ParametersDimension );
  derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  DerivativeType derivativeF = DerivativeType( ParametersDimension );
  derivativeF.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  DerivativeType derivativeM = DerivativeType( ParametersDimension );
  derivativeM.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  ti.GoToBegin();
  // First compute the sums
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
      sff += fixedValue  * fixedValue;
      smm += movingValue * movingValue;
      sfm += fixedValue  * movingValue;
      if ( this->m_SubtractMean )
        {
        sf += fixedValue;
        sm += movingValue;
        }
      this->m_NumberOfPixelsCounted++;
      }

    ++ti;
    }

  // Compute contributions to derivatives
  ti.GoToBegin();
  while(!ti.IsAtEnd())
    {

    index = ti.GetIndex();
    
    typename Superclass::InputPointType inputPoint;
    fixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

    if ( this->m_FixedImageMask && !this->m_FixedImageMask->IsInside( inputPoint ) )
      {
      ++ti;
      continue;
      }

    typename Superclass::OutputPointType transformedPoint = this->m_Transform->TransformPoint( inputPoint );

    if ( this->m_MovingImageMask && !this->m_MovingImageMask->IsInside( transformedPoint ) )
      {
      ++ti;
      continue;
      }

    if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
      {
      const RealType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );
      const RealType fixedValue     = ti.Get();

      const TransformJacobianType & jacobian =
        this->m_Transform->GetJacobian( inputPoint ); 

      // Get the gradient by NearestNeighboorInterpolation: 
      // which is equivalent to round up the point components.
      typedef typename Superclass::OutputPointType OutputPointType;
      typedef typename OutputPointType::CoordRepType CoordRepType;
      typedef ContinuousIndex<CoordRepType,MovingImageType::ImageDimension>
        MovingImageContinuousIndexType;

      MovingImageContinuousIndexType tempIndex;
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex( transformedPoint, tempIndex );

      typename MovingImageType::IndexType mappedIndex; 
      mappedIndex.CopyWithRound( tempIndex );
      
      const GradientPixelType gradient = 
        this->GetGradientImage()->GetPixel( mappedIndex );

      for(unsigned int par=0; par<ParametersDimension; par++)
        {
        RealType sumF = NumericTraits< RealType >::Zero;
        RealType sumM = NumericTraits< RealType >::Zero;
        for(unsigned int dim=0; dim<dimension; dim++)
          {
          const RealType differential = jacobian( dim, par ) * gradient[dim];
          sumF += fixedValue  * differential;
          sumM += movingValue * differential;
          if ( this->m_SubtractMean && this->m_NumberOfPixelsCounted > 0 )
            {
            sumF -= differential * sf / this->m_NumberOfPixelsCounted;
            sumM -= differential * sm / this->m_NumberOfPixelsCounted;
            }
          }
        derivativeF[par] += sumF;
        derivativeM[par] += sumM;
        }
      }

    ++ti;
    }

  if ( this->m_SubtractMean && this->m_NumberOfPixelsCounted > 0 )
    {
    sff -= ( sf * sf / this->m_NumberOfPixelsCounted );
    smm -= ( sm * sm / this->m_NumberOfPixelsCounted );
    sfm -= ( sf * sm / this->m_NumberOfPixelsCounted );
    }

  const RealType denom = -1.0 * vcl_sqrt(sff * smm );

  if( this->m_NumberOfPixelsCounted > 0 && denom != 0.0)
    {
    for(unsigned int i=0; i<ParametersDimension; i++)
      {
      derivative[i] = ( derivativeF[i] - (sfm/smm)* derivativeM[i] ) / denom;
      }
    }
  else
    {
    for(unsigned int i=0; i<ParametersDimension; i++)
      {
      derivative[i] = NumericTraits< MeasureType >::Zero;
      }
    }

  std::cout << "Derivative: " << derivative[0] << "  " << derivative[1] << "  " << derivative[2] << "  " << std::endl;  

}


// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  > void
SamplesNormalizedCorrelationImageToImageMetric<TFixedImage,TMovingImage>::ReinitializeSeed()
{
  // This method should be the same used in the ImageRandomIterator
  //vnl_sample_reseed();
Statistics::MersenneTwisterRandomVariateGenerator::GetInstance()->SetSeed();
}

// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  > void
SamplesNormalizedCorrelationImageToImageMetric<TFixedImage,TMovingImage>
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
SamplesNormalizedCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::SetNumberOfSamples(unsigned long num) const 
{
m_NumberOfSpatialSamples = num;
}



/*
 * Get both the match Measure and theDerivative Measure 
 */
template <class TFixedImage, class TMovingImage> 
void
SamplesNormalizedCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::GetValueAndDerivative(const TransformParametersType & parameters, 
                        MeasureType & value, DerivativeType  & derivative) const
{


  if( !this->GetGradientImage() )
    {
    itkExceptionMacro(<<"The gradient image is null, maybe you forgot to call Initialize()");
    }

  FixedImageConstPointer fixedImage = this->m_FixedImage;

  if( !fixedImage ) 
    {
    itkExceptionMacro( << "Fixed image has not been assigned" );
    }

  const unsigned int dimension = FixedImageType::ImageDimension;

  typedef  itk::ImageRandomConstIteratorWithIndex<FixedImageType> FixedIteratorType;

  FixedIteratorType ti( fixedImage, this->GetFixedImageRegion() );
  ti.SetNumberOfSamples(m_NumberOfSpatialSamples);

  typename FixedImageType::IndexType index;

  this->m_NumberOfPixelsCounted = 0;

  this->SetTransformParameters( parameters );

  typedef  typename NumericTraits< MeasureType >::AccumulateType AccumulateType;

  AccumulateType sff  = NumericTraits< AccumulateType >::Zero;
  AccumulateType smm  = NumericTraits< AccumulateType >::Zero;
  AccumulateType sfm = NumericTraits< AccumulateType >::Zero;
  AccumulateType sf  = NumericTraits< AccumulateType >::Zero;
  AccumulateType sm  = NumericTraits< AccumulateType >::Zero;

  const unsigned int ParametersDimension = this->GetNumberOfParameters();
  derivative = DerivativeType( ParametersDimension );
  derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  DerivativeType derivativeF = DerivativeType( ParametersDimension );
  derivativeF.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  DerivativeType derivativeM = DerivativeType( ParametersDimension );
  derivativeM.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  DerivativeType derivativeM1 = DerivativeType( ParametersDimension );
  derivativeM1.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

  ti.GoToBegin();
  // First compute the sums
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
      sff += fixedValue  * fixedValue;
      smm += movingValue * movingValue;
      sfm += fixedValue  * movingValue;
      if ( this->m_SubtractMean )
        {
        sf += fixedValue;
        sm += movingValue;
        }
      this->m_NumberOfPixelsCounted++;
      }

    ++ti;
    }


  // Compute contributions to derivatives
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
      const RealType fixedValue     = ti.Get();

      const TransformJacobianType & jacobian =
        this->m_Transform->GetJacobian( inputPoint ); 

      // Get the gradient by NearestNeighboorInterpolation: 
      // which is equivalent to round up the point components.
      typedef typename Superclass::OutputPointType OutputPointType;
      typedef typename OutputPointType::CoordRepType CoordRepType;
      typedef ContinuousIndex<CoordRepType,MovingImageType::ImageDimension>
        MovingImageContinuousIndexType;

      MovingImageContinuousIndexType tempIndex;
      this->m_MovingImage->TransformPhysicalPointToContinuousIndex( transformedPoint, tempIndex );

      typename MovingImageType::IndexType mappedIndex; 
      mappedIndex.CopyWithRound( tempIndex );
      
      const GradientPixelType gradient = 
        this->GetGradientImage()->GetPixel( mappedIndex );

      for(unsigned int par=0; par<ParametersDimension; par++)
        {
        RealType sumF = NumericTraits< RealType >::Zero;
        RealType sumM = NumericTraits< RealType >::Zero;
        for(unsigned int dim=0; dim<dimension; dim++)
          {
          const RealType differential = jacobian( dim, par ) * gradient[dim];
          sumF += fixedValue  * differential;
          sumM += movingValue * differential;
          if ( this->m_SubtractMean && this->m_NumberOfPixelsCounted > 0 )
            {
            sumF -= differential * sf / this->m_NumberOfPixelsCounted;
            sumM -= differential * sm / this->m_NumberOfPixelsCounted;
            }
          }
        derivativeF[par] += sumF;
        derivativeM[par] += sumM;
        }
      }
    ++ti;
    }

  if ( this->m_SubtractMean && this->m_NumberOfPixelsCounted > 0 )
    {
    sff -= ( sf * sf / this->m_NumberOfPixelsCounted );
    smm -= ( sm * sm / this->m_NumberOfPixelsCounted );
    sfm -= ( sf * sm / this->m_NumberOfPixelsCounted );
    }

  const RealType denom = -1.0 * vcl_sqrt(sff * smm );

  if( this->m_NumberOfPixelsCounted > 0 && denom != 0.0)
    {
    for(unsigned int i=0; i<ParametersDimension; i++)
      {
      derivative[i] = ( derivativeF[i] - (sfm/smm)* derivativeM[i] ) / denom;
      }
    value = sfm / denom;
    }
  else
    {
    for(unsigned int i=0; i<ParametersDimension; i++)
      {
      derivative[i] = NumericTraits< MeasureType >::Zero;
      }
    value = NumericTraits< MeasureType >::Zero;
    }

  
  std::cout << "Normalized correlation = " << -value << " ->  ";
//std::cout << "Normalized correlation = " << -value << " (Derivative: " << derivative[0] << "  " << derivative[1] << "  " << derivative[2] << " )  ->  ";  

}

template < class TFixedImage, class TMovingImage> 
void
SamplesNormalizedCorrelationImageToImageMetric<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "SubtractMean: " << m_SubtractMean << std::endl;
}

} // end namespace itk


#endif
