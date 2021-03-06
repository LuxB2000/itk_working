#ifndef __TRANSLATIONNORMMULTIREGISTRATION__
#define __TRANSLATIONNORMMULTIREGISTRATION__
/*
 * ===================================================================================
 *
 *       Filename:  affineMultiRegistration.h
 *
 *    Description:  Object that compute an affine registration
 *
 *        Version:  1.0
 *        Created:  06/08/14 11:28:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * ===================================================================================
 */


// registration
#include "itkMultiResolutionImageRegistrationMethod.h"

// metric
#include "itkMutualInformationImageToImageMetric.h"

// interpolator
#include "itkLinearInterpolateImageFunction.h"

// transform
#include "itkTranslationTransform.h"

// optimizer
#include "itkGradientDescentOptimizer.h"

// command iteration update
#include "translationCommandIterationUpdate.h"

// observers
#include "simpleGradientObserver.h"

template <class TMovingImage, class TFixedImage, const unsigned int Dimension = 3>
class AffineMultiRegistration{
public:

	//
	// public types
	//

	typedef TMovingImage MovingImageType;
	typedef TFixedImage  FixedImageType;
	
	// transform
  typedef itk::TranslationTransform< double, Dimension>			TransformType;

	// optimizer
  typedef itk::GradientDescentOptimizer OptimizerType;

	// metric
  typedef itk::MutualInformationImageToImageMetric< 
                                    MovingImageType, 
                                    MovingImageType >   MetricType;

	// interpolator
  typedef itk:: LinearInterpolateImageFunction< 
                                    MovingImageType,
                                    double          >    InterpolatorType;

	// registration
  typedef itk::MultiResolutionImageRegistrationMethod< 
                                    FixedImageType, 
                                    MovingImageType >   RegistrationType;

  typedef itk::MultiResolutionPyramidImageFilter<
                                    FixedImageType,
                                    FixedImageType > FixedImagePyramidType;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    MovingImageType,
                                    MovingImageType >MovingImagePyramidType;

	//
	// public functions
	//
	AffineMultiRegistration();
	~AffineMultiRegistration();

	void SetMovingImage( typename MovingImageType::Pointer moving){
		m_movingImage = moving;
	}
	void SetFixedImage( typename FixedImageType::Pointer fixed){
		m_fixedImage = fixed;
	}
	void SetPixelRatio( double ratio ){
		m_pixelRatio = ratio;
	}
	void SetNumberOfLevels( int nbrOfLevels ){
		m_nbrOfLevels = nbrOfLevels;
	}
	void SetMomentScale( double momentS ){
		m_momentScale = momentS;
	}
	void SetTranslationScale( double translationS ){
		m_translationScale = translationS;
	}
	void SetInitTranslation( float Tx, float Ty, float Tz ){
		m_initTranslationX = Tx;
		m_initTranslationY = Ty;
		m_initTranslationZ = Tz;
	}
	void SetMaxIter( int maxIter ){
		m_maxIter = maxIter;
	}
	void SetVisualization( bool v ){
		m_visualize = v;
	}

	void Update();

	OptimizerType::ParametersType GetFinalParameters(){
		return m_finalParameters;
	}
	OptimizerType::ParametersType GetFixedParameters(){
		return m_fixedParameters;
	}



	//
	// Private
	//
private:
	typename MovingImageType::Pointer m_movingImage;
	typename FixedImageType::Pointer m_fixedImage;
	double m_pixelRatio, m_momentScale, m_translationScale;
	int m_nbrOfLevels, m_maxIter;
	float m_initTranslationX, m_initTranslationY, m_initTranslationZ;

  OptimizerType::ParametersType m_finalParameters, m_fixedParameters;
	bool m_visualize;
};



#include "translationNormalisedMultiRegistration.cxx"
#endif
