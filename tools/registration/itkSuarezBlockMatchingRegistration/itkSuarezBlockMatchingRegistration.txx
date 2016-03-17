/*=========================================================================

  This work has been made to be used with the libraries Insight Segmentation
  & Registration Toolkit 2.8.1.

  Program:    itkSuarezBlockMatchingRegistration
  Module:     itkSuarezBlockMatchingRegistration.txx
  Created by: Eduardo Suarez, Rafael Nebot
              Instituto Tecnologico de Canarias - Gobierno de Canarias
              http://www.itccanarias.org/
  Date:       2006/08/07


  License:

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

// OPT: Optimize this
// DEB: Debugging info
// TOD: To do yet

// PROBAR 3D

#ifndef _itkSuarezBlockMatchingRegistration_txx
#define _itkSuarezBlockMatchingRegistration_txx

#include "itkSuarezBlockMatchingRegistration.h"
#include "itkVectorDivideImageFilter.h"
#include "itkNanWarpImageFilter.h"

#include <vnl/algo/vnl_determinant.h>
#include <vnl/vnl_trace.h>
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkVectorResampleImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkDiscreteGaussianImageFilter.h" // it has UseImageSpacingOff
//#include "itkSmoothingRecursiveGaussianImageFilter.h"

namespace itk
{

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::SuarezBlockMatchingRegistration()
{
  m_VarianceSimilarity = 2.0;
  m_VarianceRegularization = 2.0;
  m_Neighborhood.SetRadius( 1 );
  m_NumberOfLevelsNotToCompute = 0;

} // end SuarezBlockMatchingRegistration

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::VectorImagePointer
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::VectorDiscreteSmoothing( VectorImageConstPointer vectorData,
			   Real                    variance )
{
  typedef VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType>
    indexSelectionType;
  typename indexSelectionType::Pointer indexSelector =
    indexSelectionType::New();
  indexSelector->SetInput( vectorData );

  typedef DiscreteGaussianImageFilter<ScalarImageType, ScalarImageType>
    ScalarSmoothFilterType;
  typename ScalarSmoothFilterType::Pointer gaussianFilter
    = ScalarSmoothFilterType::New();

  VectorImagePointer result = VectorImageType::New();
  result->CopyInformation( vectorData );
  result->SetRegions( vectorData->GetLargestPossibleRegion() );
  result->Allocate();

  for ( unsigned index = 0; index < ImageDimension; ++index )
    {
    indexSelector->SetIndex( index );
    indexSelector->Update();
    
    gaussianFilter->SetUseImageSpacingOff();
    gaussianFilter->SetVariance( variance );
    gaussianFilter->SetInput( indexSelector->GetOutput() );
    gaussianFilter->Update();

    ScalarImagePointer component = gaussianFilter->GetOutput();

    ImageRegionIterator<VectorImageType> vectorIterator = 
      ImageRegionIterator<VectorImageType>(
	result,
	result->GetLargestPossibleRegion() );
    
    ImageRegionConstIterator<ScalarImageType> componentIterator =
      ImageRegionConstIterator<ScalarImageType>(
	component,
 	component->GetLargestPossibleRegion() );
    
    vectorIterator.GoToBegin();
    componentIterator.GoToBegin();
    
    while ( ! vectorIterator.IsAtEnd() )
      {
      vectorIterator.Value()[ index ] = componentIterator.Get();
      ++vectorIterator;
      ++componentIterator;
      } // end while

    } // end for index

  return result;

} // end VectorDiscreteSmoothing

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::ScalarImagePointer
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::LinearResampleIntoX( ScalarImageConstPointer scalarData,
		       VectorImageConstPointer deformationInX )
{
  typedef LinearInterpolateImageFunction<ScalarImageType, double>
    LinearInterpolatorType; // TOD: can't take 'double' off with gcc
  typename LinearInterpolatorType::Pointer linearInterpolator =
    LinearInterpolatorType::New();

  typedef itk::NanWarpImageFilter
    <ScalarImageType, ScalarImageType, VectorImageType> WarpType;
  typename WarpType::Pointer warp = WarpType::New();
  
  warp->SetInterpolator( linearInterpolator );
  warp->SetEdgePaddingValue(
    NumericTraits<typename WarpType::PixelType>::quiet_NaN() );
  warp->SetOutputSpacing( deformationInX->GetSpacing() );
  warp->SetOutputOrigin( deformationInX->GetOrigin() );
  warp->SetDeformationField( deformationInX );
  warp->SetInput( scalarData );
  warp->Update();
  
  return warp->GetOutput();

} // end LinearResampleIntoX (scalar)

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::VectorImagePointer
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::LinearResampleIntoX( VectorImageConstPointer vectorData,
		       VectorImageConstPointer deformationInX )
{
  VectorImagePointer result = VectorImageType::New();
  result->SetRegions( deformationInX->GetLargestPossibleRegion() );
  result->CopyInformation( deformationInX );
  result->Allocate();
  result->FillBuffer( NumericTraits<VectorPixelType>::Zero );
  
  // Filter type to split a vector in its components ///////////////////////////
  typedef VectorIndexSelectionCastImageFilter<VectorImageType, ScalarImageType>
    indexSelectionType;

  typename indexSelectionType::Pointer indexSelector =
    indexSelectionType::New();
  indexSelector->SetInput( vectorData );
  
  // Sample each component separately //////////////////////////////////////////
  for ( unsigned index = 0; index < ImageDimension; ++index )
    {
    indexSelector->SetIndex( index );
    indexSelector->Update();
    
    ScalarImagePointer component = LinearResampleIntoX(
      (ScalarImageConstPointer) indexSelector->GetOutput(), 
      deformationInX );
    
    ImageRegionIterator<VectorImageType> vectorIterator = 
      ImageRegionIterator<VectorImageType>(
	result,
	result->GetLargestPossibleRegion() );
    
    ImageRegionConstIterator<ScalarImageType> componentIterator =
      ImageRegionConstIterator<ScalarImageType>(
	component,
 	component->GetLargestPossibleRegion() );
    
    vectorIterator.GoToBegin();
    componentIterator.GoToBegin();
    
    while ( ! vectorIterator.IsAtEnd() )
      {
      vectorIterator.Value()[ index ] = componentIterator.Get();
      ++vectorIterator;
      ++componentIterator;
      } // end while

    } // end for component

  return result;

} // end LinearResampleIntoX (vector)

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::VectorImagePointer
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::ResampleToSize( VectorImageConstPointer vectorData,
		  const VectorSizeType    &size )
{
  // interpolator
  typedef VectorLinearInterpolateImageFunction< VectorImageType, double>
    LinearInterpolatorType;
  typename LinearInterpolatorType::Pointer interpolatorLinear =
    LinearInterpolatorType::New();
  
  // transformation
  typedef IdentityTransform<double, ImageDimension>
    TransformType;
  typename TransformType::Pointer transform =
    TransformType::New();

  // ivars
  Vector<double, ImageDimension> scaling;
  VectorSpacingType spacing;
  for (unsigned coordinate = 0; coordinate < ImageDimension; ++coordinate)
    {
    scaling[ coordinate ] = size[ coordinate ] /
      vectorData->GetLargestPossibleRegion().GetSize( coordinate );

    spacing[ coordinate ] =
      vectorData->GetSpacing()[ coordinate ] / scaling[ coordinate ];
    }

  // default pixel
  VectorPixelType nanVector;
  for ( int coordinate = 0; coordinate < ImageDimension; ++coordinate )
    nanVector[ coordinate ] = NumericTraits<ScalarPixelType>::quiet_NaN();

  // resampler
  typedef VectorResampleImageFilter<VectorImageType, VectorImageType>
    ResampleType;
  typename ResampleType::Pointer resampler = ResampleType::New();
  resampler->SetInterpolator( interpolatorLinear );
  resampler->SetTransform( transform );
  resampler->SetDefaultPixelValue( nanVector );
  resampler->SetOutputSpacing( spacing );
  resampler->SetOutputOrigin( vectorData->GetOrigin() );
  resampler->SetSize( size );  
  resampler->SetInput( vectorData );
  resampler->Update();
  
  return resampler->GetOutput();

} // end ResampleToSize

//------------------------------------------------------------------------------

template< class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::ScalarImagePointer
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::Structure( ScalarImageConstPointer imageData )
{
  // Compute the gradient of input image
  typename GradientImageFilter<ScalarImageType, Real, Real>::Pointer    
    gradient = GradientImageFilter<ScalarImageType, Real, Real>::New();
  gradient->SetUseImageSpacingOff();
  gradient->SetInput( imageData );
  gradient->Update();

  // Compute the operation GRAD(Image)*GRAD(Image)'
  typename UnaryFunctorImageFilter
    <CovariantVectorImageType, MatrixImageType, Mult>::Pointer internalProduct =
    UnaryFunctorImageFilter<CovariantVectorImageType, MatrixImageType, Mult>
    ::New();
  internalProduct->SetInput( gradient->GetOutput() );
  internalProduct->Update();

  // Convolve the last result with an scalar mask
  // Input image iterator
  typedef ConstNeighborhoodIterator<MatrixImageType> NeighborhoodIteratorType;

  // Output image iterator
  ScalarImagePointer output = ScalarImageType::New();
  output->SetRegions( internalProduct->GetOutput()->GetLargestPossibleRegion() );
  output->CopyInformation( internalProduct->GetOutput() );
  output->Allocate();

  ImageRegionIterator<ScalarImageType> outputIterator(
    output,
    internalProduct->GetOutput()->GetRequestedRegion());

  NeighborhoodRadiusType radius;
  radius.Fill(1);
  NeighborhoodIteratorType neighborhoodIterator(
    radius,
    internalProduct->GetOutput(), 
    internalProduct->GetOutput()->GetRequestedRegion() );
  
  // Create an array with weights corresponding to each cell of the mask
  unsigned numberOfDirections = 1;
  for ( int coordinate = 0; coordinate < ImageDimension; ++coordinate )
    numberOfDirections *= 3;

  double *weights = new double[ numberOfDirections ];
  
  // hack to compute weigths ///////////////////////////////////////////////////

  // Weight of elements: All zeros (center)= 0.20,
  //                     All (+/-) ones = 0.2/(2^D),
  //                     Others= 0.6/(3^D-2^D-1). 

  // Number of offsets wich will have ones (+1 or -1) in every coordinate
  // equals 2^ImageDimension (corners of the mask)
  unsigned nAllOnes = 1 << ImageDimension;

  double wAllZeros = 0.20 / 1;
  double wAllOnes = 0.20 / nAllOnes;
  double wMixedZerosAndOnes = 0.6 / (numberOfDirections-nAllOnes-1);

  // weights assignments to each direction
  double acum = 0.0;
  for ( unsigned direction=0; direction<numberOfDirections; ++direction )
    {
    NeighborhoodOffsetType offs = m_Neighborhood.GetOffset(direction);

    bool allZeros = true, allOnes = true;

    // guess the type of weight
    for ( unsigned coordinate=0;
	  coordinate < offs.GetOffsetDimension();
	  ++coordinate )
      {
      if (offs[coordinate]!=0)
	allZeros = false;
      else
	allOnes = false;
      }

    if (allZeros)
      weights[direction] = wAllZeros;
    else if (allOnes)
      weights[direction] = wAllOnes;
    else
      weights[direction] = wMixedZerosAndOnes;

    acum += weights[direction];
    }

  //should sum to one
  //for (unsigned i=0; i<numberOfDirections; ++i)
  //  weights[i] /= acum;
  
  for ( 
    neighborhoodIterator.GoToBegin(), outputIterator.GoToBegin();
    ! neighborhoodIterator.IsAtEnd(); 
    ++neighborhoodIterator, ++outputIterator )
    {
    MatrixType m;
    m.fill(0);

    // compute the convolution
    for ( unsigned direction=0; 
	  direction < numberOfDirections; 
	 ++direction )
      m += neighborhoodIterator.GetPixel( direction ) *
	(typename MatrixType::element_type) weights[ direction ];

    // Compute determinant of "m" & divide by its trace
    double trace = vnl_trace(m);
    if (fabs(trace)<1e-10)
      outputIterator.Set(0.0);
    else
      outputIterator.Set(vnl_determinant(m) / trace);

    }
  delete weights;
  
  return output;

} // end Structure

//------------------------------------------------------------------------------

template< class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::ScalarImagePointer 
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::SmoothWithCertainty(
  ScalarImageConstPointer signal,
  ScalarImageConstPointer certainty,
  Real                    gaussianVariance )
{
  // data * certainty //////////////////////////////////////////////////////////
  typedef MultiplyImageFilter<ScalarImageType, ScalarImageType, ScalarImageType>
    MultiplyFilterType;
  typename MultiplyFilterType::Pointer multiplyFilter =
    MultiplyFilterType::New();
  multiplyFilter->SetInput1( signal );
  multiplyFilter->SetInput2( certainty );

  // smooth( data * certainty ) ////////////////////////////////////////////////
  typedef DiscreteGaussianImageFilter<ScalarImageType, ScalarImageType>
    ScalarSmoothFilterType;
  typename ScalarSmoothFilterType::Pointer gaussianFilterDC
    = ScalarSmoothFilterType::New();
  gaussianFilterDC->SetUseImageSpacingOff();
  gaussianFilterDC->SetVariance( gaussianVariance );
  gaussianFilterDC->SetInput( multiplyFilter->GetOutput() );

  // smooth( certainty ) ///////////////////////////////////////////////////////
  typename ScalarSmoothFilterType::Pointer gaussianFilterC
    = ScalarSmoothFilterType::New();
  gaussianFilterC->SetUseImageSpacingOff();
  gaussianFilterC->SetVariance( gaussianVariance );
  gaussianFilterC->SetInput( certainty );

  // smooth( data * certainty ) / smooth( certainty ) //////////////////////////
  typedef DivideImageFilter<ScalarImageType, ScalarImageType, ScalarImageType>
    DivideFilterType;
  typename DivideFilterType::Pointer divideFilter =
    DivideFilterType::New();
  divideFilter->SetInput1( gaussianFilterDC->GetOutput() );
  divideFilter->SetInput2( gaussianFilterC->GetOutput() );

  divideFilter->Update();

  return divideFilter->GetOutput();

} // end SmoothWithCertainty (scalar)

//------------------------------------------------------------------------------

template< class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::VectorImagePointer 
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::SmoothWithCertainty(
  VectorImageConstPointer signal,
  ScalarImageConstPointer certainty,
  Real                    gaussianVariance )
{
  // data * certainty //////////////////////////////////////////////////////////
  typedef MultiplyImageFilter<VectorImageType, ScalarImageType, VectorImageType>
    MultiplyFilterType;
  typename MultiplyFilterType::Pointer multiplyFilter =
    MultiplyFilterType::New();
  multiplyFilter->SetInput1( signal );
  multiplyFilter->SetInput2( certainty );
  multiplyFilter->Update();


  // smooth( data * certainty ) ////////////////////////////////////////////////
  VectorImagePointer smoothedDataAndCertainty = VectorDiscreteSmoothing(
    multiplyFilter->GetOutput(),
    gaussianVariance );

  // smooth( certainty ) ///////////////////////////////////////////////////////
  typedef DiscreteGaussianImageFilter<ScalarImageType, ScalarImageType>
    ScalarSmoothFilterType;
  typename ScalarSmoothFilterType::Pointer gaussianFilterC 
    = ScalarSmoothFilterType::New();
  gaussianFilterC->SetUseImageSpacingOff();
  gaussianFilterC->SetVariance( gaussianVariance );
  gaussianFilterC->SetInput( certainty );

  // smooth( data * certainty ) / smooth( certainty ) //////////////////////////
  typedef VectorDivideImageFilter
    <VectorImageType, ScalarImageType, VectorImageType> VectorDivideFilterType;
  typename VectorDivideFilterType::Pointer divideFilter =
    VectorDivideFilterType::New();
  divideFilter->SetInput1( smoothedDataAndCertainty );
  divideFilter->SetInput2( gaussianFilterC->GetOutput() );
  divideFilter->Update();

  return divideFilter->GetOutput();

} // end SmoothWithCertainty (vector)


//------------------------------------------------------------------------------

template< class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::ScalarImagePointer 
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::StripNanImage( ScalarImagePointer signalWithNans )
{
  ScalarImagePointer nanImage = ScalarImageType::New();
  nanImage->CopyInformation( signalWithNans );
  nanImage->SetRegions( signalWithNans->GetLargestPossibleRegion() );
  nanImage->Allocate();
  nanImage->FillBuffer( NumericTraits<ScalarPixelType>::One );

  ImageRegionIterator<ScalarImageType> signalIterator =
    ImageRegionIterator<ScalarImageType>(
      signalWithNans,
      signalWithNans->GetLargestPossibleRegion() );

  ImageRegionIterator<ScalarImageType> nanIterator =
    ImageRegionIterator<ScalarImageType>(
      nanImage,
      nanImage->GetLargestPossibleRegion() );

  signalIterator.GoToBegin();
  nanIterator.GoToBegin();
  while ( ! nanIterator.IsAtEnd() )
    {
    if ( signalIterator.Get() != signalIterator.Get() ) // that is, if it is nan
      {
      nanIterator.Set(    NumericTraits<ScalarPixelType>::Zero );
      signalIterator.Set( NumericTraits<ScalarPixelType>::Zero );
      }
    ++nanIterator;
    ++signalIterator;
    } // end while

  return nanImage;

} // end StripNanImage (scalar)

//------------------------------------------------------------------------------

template< class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::ScalarImagePointer 
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::StripNanImage( VectorImagePointer signalWithNans )
{
  ScalarImagePointer nanImage = ScalarImageType::New();
  nanImage->CopyInformation( signalWithNans );
  nanImage->SetRegions( signalWithNans->GetLargestPossibleRegion() );
  nanImage->Allocate();
  nanImage->FillBuffer( NumericTraits<ScalarPixelType>::One );

  ImageRegionIterator<VectorImageType> signalIterator =
    ImageRegionIterator<VectorImageType>(
      signalWithNans,
      signalWithNans->GetLargestPossibleRegion() );

  ImageRegionIterator<ScalarImageType> nanIterator =
    ImageRegionIterator<ScalarImageType>(
      nanImage,
      nanImage->GetLargestPossibleRegion() );

  signalIterator.GoToBegin();
  nanIterator.GoToBegin();

  VectorPixelType zeroVector;
  VectorPixelType vect;
  for ( int coordinate = 0; coordinate < ImageDimension; ++coordinate )
    zeroVector[ coordinate ] = 0.0;

  while ( ! nanIterator.IsAtEnd() )
    {
    vect = signalIterator.Get();
    if ( vect[0] != vect[0] ) // any is a nan
      {
      nanIterator.Set( 	NumericTraits<ScalarPixelType>::Zero );
      signalIterator.Set( zeroVector );
      }
    ++nanIterator;
    ++signalIterator;
    } // end while

  return nanImage;

} // end StripNanImage (vector)

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::ImageRegionType 
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::RegionFromDisplacementIndex(
  unsigned              displacementIndex,
  const ImageRegionType &baseRegion )
{
  ImageRegionType region;
  
  NeighborhoodOffsetType offset = m_Neighborhood.GetOffset( displacementIndex );
  
  for ( unsigned d = 0; d < ImageDimension; ++d )
    {
    // When defining the region, dismiss two cells on each
    // dimension, one to the "left" and one to the "right"
    region.SetSize( d, baseRegion.GetSize( d ) - 2 );
    region.SetIndex( d, baseRegion.GetIndex( d ) + 1 + offset[ d ] );
    } // end for
  
  return region;

} // end regionFromDirectionIndex

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::VectorType 
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::DisplacementFromDisplacementIndex( unsigned displacementIndex )
{
  NeighborhoodOffsetType offset = m_Neighborhood.GetOffset( displacementIndex );
  
  Vector<int, ImageDimension> displacement;
  for ( unsigned coordinate = 0; coordinate < ImageDimension; ++coordinate )
    displacement[ coordinate ] = offset[ coordinate ];

  return displacement;

} // end DisplacementFromDisplacementIndex

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::ScalarImagePointer
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::SimilaritySSD( unsigned                displacementIndex,
		 ScalarImageConstPointer fixedInX,
		 ScalarImageConstPointer movingInX,
		 Real                    varianceSimilarity )
{
  
  // At this point it is assumed than both Origin, Spacing and Size of
  // fixedInX and movingInX are the same.

  // WARNING: Fixed means in the reference system for computation (Image A)
  //          It does not mean it is the "fixed image". Meanings are
  //          opposite for intuitive image movement in registration
  //          and moving in computations. It is the same for Moving (B).

  unsigned zeroDisplacementIndex = m_Neighborhood.GetCenterNeighborhoodIndex();

  ImageRegionType fixedRegion = RegionFromDisplacementIndex(
    zeroDisplacementIndex,
    fixedInX->GetLargestPossibleRegion() );

  ImageRegionType movingRegion = RegionFromDisplacementIndex(
    displacementIndex,
    movingInX->GetLargestPossibleRegion() );

  ImageRegionConstIterator<ScalarImageType> fixedInXIterator =
    ImageRegionConstIterator<ScalarImageType>( fixedInX, fixedRegion );

  ImageRegionConstIterator<ScalarImageType> movingInXIterator =
    ImageRegionConstIterator<ScalarImageType>( movingInX, movingRegion );

  ScalarImagePointer ssdInX = ScalarImageType::New();
  ssdInX->SetRegions( fixedInX->GetLargestPossibleRegion() );
  ssdInX->CopyInformation( fixedInX );
  ssdInX->Allocate();
  ssdInX->FillBuffer( NumericTraits<ScalarPixelType>::quiet_NaN() );

  ImageRegionIterator<ScalarImageType> ssdInXIterator =
    ImageRegionIterator<ScalarImageType>( ssdInX, fixedRegion );

  Real similarity;
  fixedInXIterator.GoToBegin();
  movingInXIterator.GoToBegin();
  ssdInXIterator.GoToBegin();
  while ( ! ssdInXIterator.IsAtEnd() )
    {
    similarity = fixedInXIterator.Get() - movingInXIterator.Get();
    ssdInXIterator.Set( - (similarity * similarity) );

    ++fixedInXIterator;
    ++movingInXIterator;
    ++ssdInXIterator;
    } // end while

  ScalarImagePointer certaintySsdInX = StripNanImage( ssdInX );

  // TOD: smoothing with recursiveSmoothingGaussian doesn't work well
  ScalarImagePointer result = SmoothWithCertainty(
    (ScalarImageConstPointer) ssdInX,
    (ScalarImageConstPointer) certaintySsdInX,
    varianceSimilarity );

  return result;

} // end SimilaritySSD

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::ScalarImagePointer
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::SimilarityMI( unsigned                displacementIndex,
		ScalarImageConstPointer fixedInX,
		ScalarImageConstPointer movingInX,
		NanNormalizer           &normalizerFixed,
		NanNormalizer           &normalizerMoving,
		Real                    varianceSimilarity )
{
  unsigned zeroDisplacementIndex = m_Neighborhood.GetCenterNeighborhoodIndex();

  ImageRegionType fixedRegion = RegionFromDisplacementIndex(
    zeroDisplacementIndex,
    fixedInX->GetLargestPossibleRegion() );

  ImageRegionType movingRegion = RegionFromDisplacementIndex(
    displacementIndex,
    fixedInX->GetLargestPossibleRegion() );

  vnl_matrix<float> jointHistogram;
  jointHistogram.set_size(
    normalizerFixed.getNBins(),
    normalizerMoving.getNBins() );
  jointHistogram.fill( 0.0 );

  // fixed histogram makes no difference

  vnl_vector<float> movingHistogram;
  movingHistogram.set_size( normalizerMoving.getNBins() );
  movingHistogram.fill( 0.0 );

  ImageRegionConstIterator<ScalarImageType> fixedInXIterator =
    ImageRegionConstIterator<ScalarImageType>( fixedInX, fixedRegion );

  ImageRegionConstIterator<ScalarImageType> movingInXIterator =
    ImageRegionConstIterator<ScalarImageType>( movingInX, movingRegion );

  int fixedBin, movingBin;
  fixedInXIterator.GoToBegin();
  movingInXIterator.GoToBegin();
  while ( ! fixedInXIterator.IsAtEnd() )
    {
    fixedBin = normalizerFixed.getBin( fixedInXIterator.Get() );
    movingBin = normalizerMoving.getBin( movingInXIterator.Get() );
    if ( fixedBin != -1 && movingBin != -1 ) // -1 for nans
      {
      ++jointHistogram[ fixedBin ][ movingBin ];
      ++movingHistogram[ movingBin ];
      }
    ++fixedInXIterator;
    ++movingInXIterator;
    } // end while

  // histograms don't need to be normalized because they will divided
  // each other later

  ScalarImagePointer miInX = ScalarImageType::New();
  miInX->SetRegions( fixedInX->GetLargestPossibleRegion() );
  miInX->CopyInformation( fixedInX );
  miInX->Allocate();
  miInX->FillBuffer( NumericTraits<ScalarPixelType>::Zero );

  ScalarImagePointer certaintyMiInX = ScalarImageType::New();
  certaintyMiInX->CopyInformation( fixedInX );
  certaintyMiInX->SetRegions( fixedInX->GetLargestPossibleRegion() );
  certaintyMiInX->Allocate();
  certaintyMiInX->FillBuffer( NumericTraits<ScalarPixelType>::Zero );
  
  ImageRegionIterator<ScalarImageType> miInXIterator =
    ImageRegionIterator<ScalarImageType>( miInX, fixedRegion );

  ImageRegionIterator<ScalarImageType> certaintyMiInXIterator =
    ImageRegionIterator<ScalarImageType>( certaintyMiInX, fixedRegion );

  Real movingProbabilityDensity, jointProbabilityDensity;

  fixedInXIterator.GoToBegin();
  movingInXIterator.GoToBegin();
  miInXIterator.GoToBegin();
  certaintyMiInXIterator.GoToBegin();
  while ( ! fixedInXIterator.IsAtEnd() )
    {
    fixedBin = normalizerFixed.getBin(fixedInXIterator.Get());
    movingBin = normalizerMoving.getBin(movingInXIterator.Get());
    if ( fixedBin != -1 && movingBin != -1 )
      {
      movingProbabilityDensity = movingHistogram[ movingBin ];
      jointProbabilityDensity  = jointHistogram[ fixedBin ][ movingBin ];

      if ( jointProbabilityDensity > 0 )
	{
	miInXIterator.Set (
	  log( jointProbabilityDensity / movingProbabilityDensity) );
	certaintyMiInXIterator.Set ( 1.0 );
	}
      else if ( movingProbabilityDensity > 0 )
	{
	miInXIterator.Set (
	  NumericTraits<ScalarPixelType>::NonpositiveMin() );
	certaintyMiInXIterator.Set ( 1.0 );
	}
      } //end if

    ++fixedInXIterator;
    ++movingInXIterator;
    ++miInXIterator;
    ++certaintyMiInXIterator;
    } // end while

  // TOD: smoothing with recursiveSmoothingGaussian doesn't work well
  ScalarImagePointer result = SmoothWithCertainty(
    (ScalarImageConstPointer) miInX,
    (ScalarImageConstPointer) certaintyMiInX,
    varianceSimilarity );

  return result;

} // end SimilarityMI

//------------------------------------------------------------------------------

template <class TInputImage, class TOutputImage>
typename SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::VectorImagePointer
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::UpdateDeformation( ScalarImageConstPointer imageAInX,
		     ScalarImageConstPointer imageBInY,
		     VectorImagePointer      initialDeformationInX )
{

  // this will be used to create new images
  ImageRegionType bufferedRegionX( imageAInX->GetBufferedRegion() );

  // this is totally heuristic (ad hoc)
  const unsigned numberOfBins =
    static_cast<unsigned>(
      std::max( ceil( std::sqrt(
			static_cast<double>( 
			  bufferedRegionX.GetNumberOfPixels() )
			/128. ) ),
		16.0 ) );

  // A in the reference system X is the image that is deformed with
  // (that is, follows) the deformation field. So that, we work on the
  // other side and move image B in the reference system Y to the
  // reference system X, just by sampling B in the positions of the
  // deformation field defined over X.
  ScalarImagePointer imageBInX = LinearResampleIntoX(
    (ScalarImageConstPointer) imageBInY,
    (VectorImageConstPointer) initialDeformationInX );

  // Similarity image type. Set to "-infinite"
  ScalarImagePointer similarityUntilLastTestDisplacement =
    ScalarImageType::New();
  similarityUntilLastTestDisplacement->CopyInformation( imageAInX );
  similarityUntilLastTestDisplacement->SetBufferedRegion( bufferedRegionX );
  similarityUntilLastTestDisplacement->Allocate();
  similarityUntilLastTestDisplacement->FillBuffer(
    NumericTraits<ScalarPixelType>::NonpositiveMin() );
  // OPT: This last line could be avoided

  Natural numberOfTestDisplacements = 1;
  for ( int coordinate = 0; coordinate < ImageDimension; ++coordinate )
    numberOfTestDisplacements *= 3;

  NaturalScalarImagePointer lastOptimumDisplacementIndex =
    NaturalScalarImageType::New();
  lastOptimumDisplacementIndex->CopyInformation( imageAInX );
  lastOptimumDisplacementIndex->SetBufferedRegion( bufferedRegionX );
  lastOptimumDisplacementIndex->Allocate();
  lastOptimumDisplacementIndex->FillBuffer( numberOfTestDisplacements + 1 );
  
  typedef itk::MinimumMaximumImageCalculator<ScalarImageType> 
    MinMaxCalculatorType;
  typename MinMaxCalculatorType::Pointer minmaxCalculator =
    MinMaxCalculatorType::New(); //it already ignores nans

  // normalizers are used for computing histograms
  NanNormalizer normA, normB;
  minmaxCalculator->SetImage( imageAInX );
  minmaxCalculator->Compute();
  normA.setUp(
    minmaxCalculator->GetMinimum(),
    minmaxCalculator->GetMaximum(),
    numberOfBins );

  minmaxCalculator->SetImage( imageBInX );
  minmaxCalculator->Compute();
  normB.setUp(
    minmaxCalculator->GetMinimum(),
    minmaxCalculator->GetMaximum(),
    numberOfBins );

  std::cout << "/";
  std::cout.flush();
  unsigned zeroDisplacementIndex = m_Neighborhood.GetCenterNeighborhoodIndex();

  for ( unsigned testDisplacementPseudoindex=0;
	testDisplacementPseudoindex < numberOfTestDisplacements;
	++testDisplacementPseudoindex )
    {

    // permutation to start with zero displacement index ///////////////////////
    unsigned testDisplacementIndex;
    if ( testDisplacementPseudoindex == 0 )
      testDisplacementIndex = zeroDisplacementIndex;
    else if ( testDisplacementPseudoindex == zeroDisplacementIndex )
      testDisplacementIndex = 0;
    else
      testDisplacementIndex = testDisplacementPseudoindex;

    // computation of similarity ///////////////////////////////////////////////
//     ScalarImagePointer similarityOfThisTestDisplacement = SimilaritySSD(
//       testDisplacementIndex,
//       (ScalarImageConstPointer) imageAInX,
//       (ScalarImageConstPointer) imageBInX,
//       varianceSimilarity );

    ScalarImagePointer similarityOfThisTestDisplacement = SimilarityMI(
      testDisplacementIndex,
      (ScalarImageConstPointer) imageAInX,
      (ScalarImageConstPointer) imageBInX,
      normA,
      normB,
      m_VarianceSimilarity );

    // update for this level ///////////////////////////////////////////////////
    ImageRegionConstIterator<ScalarImageType>
      similarityOfThisTestDisplacementIterator =
      ImageRegionConstIterator<ScalarImageType>(
	similarityOfThisTestDisplacement,
	similarityOfThisTestDisplacement->GetLargestPossibleRegion() );

    ImageRegionIterator<ScalarImageType>
      similarityUntilLastTestDisplacementIterator =
      ImageRegionIterator<ScalarImageType>(
	similarityUntilLastTestDisplacement,
	similarityUntilLastTestDisplacement->GetLargestPossibleRegion() );

    ImageRegionIterator<NaturalScalarImageType>
      lastOptimumDisplacementIndexIterator =
      ImageRegionIterator<NaturalScalarImageType>(
	lastOptimumDisplacementIndex,
	lastOptimumDisplacementIndex->GetLargestPossibleRegion() );

    ScalarPixelType thisSimilarity;
    similarityUntilLastTestDisplacementIterator.GoToBegin();
    similarityOfThisTestDisplacementIterator.GoToBegin();
    lastOptimumDisplacementIndexIterator.GoToBegin();
    while ( ! similarityUntilLastTestDisplacementIterator.IsAtEnd() )
      {
      thisSimilarity = similarityOfThisTestDisplacementIterator.Get();

      if ( thisSimilarity > similarityUntilLastTestDisplacementIterator.Get() )
	{
	similarityUntilLastTestDisplacementIterator.Set( thisSimilarity );
	lastOptimumDisplacementIndexIterator.Set( testDisplacementIndex );
	}		

      ++similarityUntilLastTestDisplacementIterator;
      ++similarityOfThisTestDisplacementIterator;
      ++lastOptimumDisplacementIndexIterator;
      } //end while iterators

    std::cout << ".";
    std::cout.flush();

    } // end for test indices

  // turn index image into vector image ////////////////////////////////////////
  VectorImagePointer discreteDisplacementInX = VectorImageType::New();
  discreteDisplacementInX->CopyInformation( imageAInX );
  discreteDisplacementInX->SetBufferedRegion( bufferedRegionX );
  discreteDisplacementInX->Allocate();

  ImageRegionIterator<VectorImageType> discreteDisplacementInXIterator =
    ImageRegionIterator<VectorImageType>(
      discreteDisplacementInX,
      discreteDisplacementInX->GetLargestPossibleRegion() );

  ImageRegionConstIterator<NaturalScalarImageType> displacementIndexIterator =
    ImageRegionConstIterator<NaturalScalarImageType>(
      lastOptimumDisplacementIndex,
      lastOptimumDisplacementIndex->GetLargestPossibleRegion() );
  
  VectorType vect;    
  discreteDisplacementInXIterator.GoToBegin();
  displacementIndexIterator.GoToBegin();
  VectorSpacingType spacingFactor =
    discreteDisplacementInX->GetSpacing();

  while ( ! discreteDisplacementInXIterator.IsAtEnd() )
    {
    vect = DisplacementFromDisplacementIndex( displacementIndexIterator.Get() );

    for ( int coordinate = 0; coordinate < ImageDimension; ++coordinate )
      vect[ coordinate ] *= spacingFactor[ coordinate ];

    discreteDisplacementInXIterator.Set( vect );

    ++discreteDisplacementInXIterator;
    ++displacementIndexIterator;
    } // end while

  // Normalized convolution of discrete displacements //////////////////////////
  ScalarImagePointer certaintyInX = StripNanImage( imageBInX );
  
  ScalarImagePointer localStructureOfAInX = Structure (
    (ScalarImageConstPointer) imageAInX );

    typedef MultiplyImageFilter
      <ScalarImageType, ScalarImageType, ScalarImageType> MultiplyFilterType;
    typename MultiplyFilterType::Pointer multiplyFilter =
      MultiplyFilterType::New();
    multiplyFilter->SetInput1( certaintyInX );
    multiplyFilter->SetInput2( localStructureOfAInX );
    multiplyFilter->Update();

  VectorImagePointer incrementalDeformationInX = 
    SmoothWithCertainty(
      (VectorImageConstPointer) discreteDisplacementInX,
//      (ScalarImageConstPointer) certaintyInX ,
      (ScalarImageConstPointer) multiplyFilter->GetOutput() ,
      m_VarianceRegularization );

  // Composition of incremental and initial deformation ////////////////////////
  ImageRegionIterator<VectorImageType> initialDeformationInXIterator = 
    ImageRegionIterator<VectorImageType>(
      initialDeformationInX,
      initialDeformationInX->GetLargestPossibleRegion() );

  VectorPointType position;
  initialDeformationInXIterator.GoToBegin();
  while ( ! initialDeformationInXIterator.IsAtEnd() )
    {
    vect = initialDeformationInXIterator.Get();
    initialDeformationInX->TransformIndexToPhysicalPoint(
      initialDeformationInXIterator.GetIndex(),
      position );

    for ( int coordinate = 0; coordinate < ImageDimension; ++coordinate )
      vect[ coordinate ] += position[ coordinate ];

    initialDeformationInXIterator.Set( vect );
  
    ++initialDeformationInXIterator;
    } // end while
  
  VectorImagePointer dirtyFinalDeformationInX = LinearResampleIntoX(
    (VectorImageConstPointer) initialDeformationInX,
    (VectorImageConstPointer) incrementalDeformationInX );

  ImageRegionIterator<VectorImageType> dirtyFinalDeformationInXIterator = 
    ImageRegionIterator<VectorImageType>(
      dirtyFinalDeformationInX,
      dirtyFinalDeformationInX->GetLargestPossibleRegion() );
  
  dirtyFinalDeformationInXIterator.GoToBegin();
  while ( ! dirtyFinalDeformationInXIterator.IsAtEnd() )
    {
    vect = dirtyFinalDeformationInXIterator.Get();
    dirtyFinalDeformationInX->TransformIndexToPhysicalPoint(
      dirtyFinalDeformationInXIterator.GetIndex(),
      position );

    for ( int coordinate = 0; coordinate < ImageDimension; ++coordinate )
      vect[ coordinate ] -= position[ coordinate ];

    dirtyFinalDeformationInXIterator.Set( vect );
    ++dirtyFinalDeformationInXIterator;
    } // end while

  certaintyInX = StripNanImage( dirtyFinalDeformationInX );
  VectorImagePointer finalDeformationInX = 
    SmoothWithCertainty(
      (VectorImageConstPointer) dirtyFinalDeformationInX,
      (ScalarImageConstPointer) certaintyInX,
      m_VarianceRegularization );
  
  std::cout << "\\";
  std::cout.flush();

  return finalDeformationInX;
  
} // end UpdateDeformation

//------------------------------------------------------------------------------

template< class TInputImage, class TOutputImage>
void
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::GenerateData( void )
{
  std::cout << "*** Start ***" << std::endl;
  std::cout.flush();

  typename TInputImage::ConstPointer baseImageAInX = this->GetInput(0);
  typename TInputImage::ConstPointer baseImageBInY = this->GetInput(1);
  // TOD: Check that both images are initially of the same size

  typedef MultiResolutionPyramidImageFilter<InputImageType,ScalarImageType>
    PyramidType;

  const int numberOfPyramidalLevels = 9;
  const int minimumSize = 12;

  typename CommandObserver::Pointer observer = CommandObserver::New();

  typename PyramidType::Pointer pyramidAInX = PyramidType::New();
  pyramidAInX->SetInput( baseImageAInX );
  pyramidAInX->SetNumberOfLevels( numberOfPyramidalLevels ); 
  pyramidAInX->AddObserver( ProgressEvent(), observer );

  typename PyramidType::Pointer pyramidBInY = PyramidType::New();
  pyramidBInY->SetInput( baseImageBInY );
  pyramidBInY->SetNumberOfLevels( numberOfPyramidalLevels );
  pyramidBInY->AddObserver( ProgressEvent(), observer );

  // Find initial level
  int initialLevel = 0 ;
  InputSizeType size;
  do
    {
    pyramidAInX->GetOutput( initialLevel )->Update();

    size = pyramidAInX->GetOutput( initialLevel )->
      GetLargestPossibleRegion().GetSize();
    
    if ( size[ 0 ] < minimumSize ) ++initialLevel;
    }
  while ( size[ 0 ] < minimumSize );

  std::cout << std::endl

	    << "Number of pyramidal levels = " << numberOfPyramidalLevels
	    << std::endl

	    << "Initial level = " << initialLevel
	    << std::endl

	    << "Final level = " << 
    ( numberOfPyramidalLevels - m_NumberOfLevelsNotToCompute )
	    << std::endl;

  std::cout.flush();
  
  // Declare & initialize deformationInX
  VectorImagePointer deformationInX = VectorImageType::New();
  deformationInX->CopyInformation( pyramidAInX->GetOutput( initialLevel ) );
  deformationInX->SetRegions(
    pyramidAInX->GetOutput( initialLevel )->GetLargestPossibleRegion() );
  deformationInX->Allocate();
  deformationInX->FillBuffer( NumericTraits<VectorPixelType>::Zero );

  // Main loop - Iterate through each level
  VectorImagePointer auxiliarDeformationInX;
  unsigned level = initialLevel;
  for ( ;
	level < numberOfPyramidalLevels - m_NumberOfLevelsNotToCompute;
	++level )
    {
    std::cout << "Level " << level << ": ";
    std::cout.flush();

    pyramidAInX->GetOutput( level )->Update();
    std::cout << "+";
    std::cout.flush();

    pyramidBInY->GetOutput( level )->Update();
    std::cout << "+";
    std::cout.flush();


    // resampling is not necessary for the first level
    if ( level == initialLevel )
      {
      auxiliarDeformationInX = deformationInX;
      }
    else
      {
      auxiliarDeformationInX = ResampleToSize(
	(VectorImageConstPointer) deformationInX,
	( pyramidAInX->GetOutput( level )->
	  GetLargestPossibleRegion() ).GetSize() );
      }
    std::cout << "*";
    std::cout.flush();

    deformationInX = UpdateDeformation(
      pyramidAInX->GetOutput( level ),
      pyramidBInY->GetOutput( level ),
      auxiliarDeformationInX );

    std::cout << std::endl;
    std::cout.flush();

    } // end for levels
  
    if ( level == numberOfPyramidalLevels )
      {
      auxiliarDeformationInX = deformationInX;
      }
    else
      {
      auxiliarDeformationInX = ResampleToSize(
	(VectorImageConstPointer) deformationInX,
	( pyramidAInX->GetOutput( numberOfPyramidalLevels -1 )->
	  GetLargestPossibleRegion() ).GetSize() );
      }

  std::cout << "*** End ***" << std::endl;
  std::cout.flush();

  this->GraftOutput( auxiliarDeformationInX );

} // end GenerateData

//------------------------------------------------------------------------------

template< class TInputImage, class TOutputImage>
void
SuarezBlockMatchingRegistration<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

//------------------------------------------------------------------------------

} // end namespace itk

#endif



