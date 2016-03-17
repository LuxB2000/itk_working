/*=========================================================================

  This work has been made to be used against the libraries Insight Segmentation
  & Registration Toolkit 2.8.1.

  Program:    itkSuarezBlockMatchingRegistration
  Module:     itkSuarezBlockMatchingRegistration.h
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

#ifndef __itkSuarezBlockMatchingRegistration_h
#define __itkSuarezBlockMatchingRegistration_h

#include "itkImageToImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkNeighborhood.h"

#include "itkCommand.h"

namespace itk
{

/** \class SuarezBlockMatchingRegistration
 * \brief Implements a blocking matching deformable registration
 * between two images
 */

template <class TInputImage, class TOutputImage>
class ITK_EXPORT SuarezBlockMatchingRegistration :
    public ImageToImageFilter< TInputImage, TOutputImage > 
{
public:
  /** Standard "Self" & Superclass typedef. */
  typedef SuarezBlockMatchingRegistration Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SuarezBlockMatchingRegistration, ImageToImageFilter);
  
  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  typedef typename TOutputImage::PixelType         OutputPixelType;
  typedef typename TOutputImage::InternalPixelType OutputInternalPixelType;
  typedef typename  TInputImage::PixelType         InputPixelType;
  typedef typename  TInputImage::InternalPixelType InputInternalPixelType;
  typedef typename  TInputImage::SizeType          InputSizeType;

  /** Extract some information from the image types.  Dimensionality
   * of the two images is assumed to be the same. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  
  /** Image typedef support. */
  typedef TInputImage  RegisteredImageType;
  typedef TInputImage  InputImageType;
  typedef TOutputImage DisplacementImageType;
  typedef typename RegisteredImageType::Pointer RegisteredImagePointer;
  
  /** Superclass typedefs. */
  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  /** Set the moving image */
  void SetMovingImage ( RegisteredImageType* MovingImage ) {
    this->SetNthInput ( 1, MovingImage );
  }
  void SetFixedImage ( RegisteredImageType* FixedImage ) {
    this->SetInput ( FixedImage );
  }
  itkSetMacro( VarianceRegularization, double );
  itkGetConstReferenceMacro( VarianceRegularization, double );
  itkSetMacro( VarianceSimilarity, double );
  itkGetConstReferenceMacro( VarianceSimilarity, double );
  itkSetMacro( NumberOfLevelsNotToCompute, unsigned );
  itkGetConstReferenceMacro( NumberOfLevelsNotToCompute, unsigned );

protected:
  SuarezBlockMatchingRegistration();
  virtual ~SuarezBlockMatchingRegistration() {};

  void PrintSelf(std::ostream& os, Indent indent) const;
  void GenerateData();

private:

  // Types declaration /////////////////////////////////////////////////////////
  typedef float Real;
  typedef unsigned Natural;

  typedef Image<Real, TInputImage::ImageDimension> ScalarImageType;
  typedef typename ScalarImageType::Pointer ScalarImagePointer;
  typedef typename ScalarImageType::ConstPointer ScalarImageConstPointer;
  typedef typename ScalarImageType::PixelType ScalarPixelType;
  typedef typename ScalarImageType::RegionType ImageRegionType;

  typedef Image<Natural, TInputImage::ImageDimension> NaturalScalarImageType;
  typedef typename NaturalScalarImageType::Pointer NaturalScalarImagePointer;

  typedef Vector<Real, TInputImage::ImageDimension> VectorType;
  typedef Image<VectorType, TInputImage::ImageDimension> VectorImageType;
  typedef typename VectorImageType::Pointer VectorImagePointer;
  typedef typename VectorImageType::ConstPointer VectorImageConstPointer;
  typedef typename VectorImageType::PixelType VectorPixelType;
  typedef typename VectorImageType::SizeType VectorSizeType;
  typedef typename VectorImageType::PointType VectorPointType;
  typedef typename VectorImageType::SpacingType VectorSpacingType;

  typedef vnl_matrix_fixed
    <Real, TInputImage::ImageDimension, TInputImage::ImageDimension> MatrixType;
  typedef Image<MatrixType, TInputImage::ImageDimension> MatrixImageType;

  typedef typename GradientImageFilter<TInputImage, Real, Real>::OutputPixelType
    CovariantVectorType;
  typedef typename GradientImageFilter<TInputImage, Real, Real>::OutputImageType
    CovariantVectorImageType;

  typedef Neighborhood<InputPixelType, ImageDimension> NeighborhoodType;
  typedef typename NeighborhoodType::OffsetType NeighborhoodOffsetType;
  typedef typename NeighborhoodType::RadiusType NeighborhoodRadiusType;

  // Methods  //////////////////////////////////////////////////////////////////

  SuarezBlockMatchingRegistration(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  class NanNormalizer; // helper class used in 
                       // UpdateDeformation and SimilarityMI

  ImageRegionType RegionFromDisplacementIndex(
    unsigned              displacementIndex,
    const ImageRegionType &baseRegion );

  VectorType DisplacementFromDisplacementIndex(
    unsigned displacementIndex );

  VectorImagePointer VectorDiscreteSmoothing(
    VectorImageConstPointer vectorData,
    Real                    variance );

  VectorImagePointer ResampleToSize(
    VectorImageConstPointer vectorData,
    const VectorSizeType    &size );

  ScalarImagePointer LinearResampleIntoX(
    ScalarImageConstPointer scalarData,
    VectorImageConstPointer deformationInX );

  VectorImagePointer LinearResampleIntoX(
    VectorImageConstPointer vectorData,
    VectorImageConstPointer deformationInX );

  ScalarImagePointer SmoothWithCertainty(
    ScalarImageConstPointer signal,
    ScalarImageConstPointer certainty,
    Real                    gaussianVariance );

  VectorImagePointer SmoothWithCertainty(
    VectorImageConstPointer signal,
    ScalarImageConstPointer certainty,
    Real                    gaussianVariance );

  ScalarImagePointer StripNanImage(
    ScalarImagePointer signal_with_nans );

  ScalarImagePointer StripNanImage(
    VectorImagePointer signal_with_nans );

  ScalarImagePointer Structure(
    ScalarImageConstPointer imageData );

  ScalarImagePointer SimilaritySSD(
    unsigned                displacementIndex,
    ScalarImageConstPointer fixedInX,
    ScalarImageConstPointer movingInX,
    Real                    varianceSimilarity );

  ScalarImagePointer SimilarityMI(
    unsigned                displacementIndex,
    ScalarImageConstPointer fixedInX,
    ScalarImageConstPointer movingInX,
    NanNormalizer           &normalizerFixed,
    NanNormalizer           &normalizerMoving,
    Real                    varianceSimilarity );

  VectorImagePointer UpdateDeformation(
    ScalarImageConstPointer imageAInX,
    ScalarImageConstPointer imageBInY,
    VectorImagePointer      initialDeformationInX );

  NaturalScalarImagePointer normalize(
    ScalarImageConstPointer x,
    Natural bins,
    ScalarImageConstPointer cx );

  // Member variables  /////////////////////////////////////////////////////////

  double m_VarianceRegularization;
  double m_VarianceSimilarity;
  unsigned m_NumberOfLevelsNotToCompute;
  NeighborhoodType m_Neighborhood;

  // Helper Classes ////////////////////////////////////////////////////////////
  class NanNormalizer
  {
    Real m_min, m_max;
    unsigned m_binsImage;
    Real m_binWidth;
  public:
    void setUp( Real min, Real max, unsigned numberOfBins )
    {
      m_min = min;
      m_max = max;
      m_binsImage = numberOfBins;
      m_binWidth = ( m_max - m_min ) / m_binsImage;
    }
    const inline unsigned getNBins() { return m_binsImage; }
    const inline int getBin( const Real &value )
    {
      // bins index: 0..m_binsIndex-1 for scalar range
      //                           -1 for nans
      return value != value ? -1 :
	( value <= m_min ? 0 :
	  ( value >= m_max ? m_binsImage-1 : 
	    static_cast<unsigned>( ( value - m_min ) / m_binWidth ) ) );
    }
  }; // end class NanNormalizer
  
  class Mult
  {
  public:
    Mult() {}
    ~Mult() {}
    inline MatrixType operator()( const CovariantVectorType & A )
    {
      // Multiply A (column vector) by its transpose (row vector)
      MatrixType m;
      for ( unsigned i = 0; i < A.GetCovariantVectorDimension(); ++i )
	for ( unsigned j = i; j < A.GetCovariantVectorDimension(); ++j )
	  m[ i ][ j ] = m[ j ][ i ] = A[ i ] * A[ j ];
      return m;
    }
  }; // end class Mult
  
  class CommandObserver : public Command
  {
  public:
    typedef  CommandObserver          Self;
    typedef  Command             Superclass;
    typedef  SmartPointer<Self>  Pointer;
    itkNewMacro( Self );
  protected:
    CommandObserver() {};
  public:
    typedef   itk::ProcessObject      ProcessType;
    typedef   const ProcessType   *   ProcessPointer;

    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      ProcessPointer filter = dynamic_cast< ProcessPointer >( object );
      if( typeid( event ) == typeid( itk::ProgressEvent ) )
	{
	std::cout << "p";
	std::cout.flush();
	}
    }
  };

}; // end class SuarezBlockMatchingRegistration

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSuarezBlockMatchingRegistration.txx"
#endif

#endif

