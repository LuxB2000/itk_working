#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

/*
 * original authors: Jonathan Orban de Xirvy and Guillaume Janssens, UCL
 * modified by Jerome Plumat, UoA
 * modifications: optimizers scales have been modified: the rotation scales
 * are based on the spacings along dimensions and translation are constant
 * EXCEPT the translation along Z that is bigger due to GP application
 */

#include <stdlib.h>
#include <stdio.h>
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkShrinkImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkImage.h"
//#include "itkOrientedImage.h"
#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkVersorRigid3DTransformOptimizer.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMetaImageIO.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkResampleImageFilter.h"
#include "itkCommand.h"
//#include "itkMattesMutualInformationImageToImageMetricNew.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkNormalizeImageFilter.h"
#include "itkCastImageFilter.h"
#include "vnl/vnl_math.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkTransformFileWriter.h"
#include "itkHistogramMatchingImageFilter.h"
// Class commanditerationupdate
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:

template <class TOptimizer>
class CommandIterationUpdate : public itk::Command {
public:
  typedef CommandIterationUpdate Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef const TOptimizer * OptimizerPointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() { }
  ~CommandIterationUpdate() { }
    // the non-const execute mehtod just calls the const one
  void Execute( itk::Object * caller, const itk::EventObject & event ){
    Execute( (const itk::Object *)caller, event );}
   // const execute method
  void Execute( const itk::Object * object, const itk::EventObject & event ){
    // get the optimizer that triggered the event
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
    // was there really an iteration Event?
    if ( ! itk::IterationEvent().CheckEvent( &event ) ) return;
    // print iteration and value to screen

    std::cout << optimizer->GetCurrentIteration();
		std::cout << " = " << (double)(optimizer->GetValue()) << std::endl;
}
};

// Main
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
	std::cout<<" -- "<<argv[0]<<" -- "<<std::endl;
  if( (argc<6) || (argc>11) ){
    std::cerr << "Wrong number of Parameters (" << argc << ")" << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " FixedImageDir MovingImageDir OutputImageDir Number_of_Levels Min_Number_of_Iterations\n" 
			<< " (outputTransformFile) (outputInverseTransform) (inputTransform=[x,y,z]; default [0 0 0])" << std::endl << std::endl;
    std::cerr << "Provided parameters :" << std::endl;
	for (int i = 1; i < argc; i++){
		std::cout << "ARG[" << i << "] : " << argv[i] << std::endl;}
    return EXIT_FAILURE;
  }
	// ARGC
	const unsigned int 
			argcFixedImage = 1, argcMovingImage = 2, argcOutputImage = 3,
			argcNbrOfLevels = 4, argcNbrOfIterations = 5,
			argcOutputTransf = 6, argcOutputInvTransf = 7, 
			argcInputTransfomX = 8, argcInputTransfomY = 9, argcInputTransfomZ = 10;

  std::ofstream OutputTransform, OutputInvTransform;
  const char* OutputTransformFileName;
	const char* OutputInvTransformFileName;
	float inputTransformX = 0, inputTransformY = 0, inputTransformZ = 0;
  if(argc >= (argcOutputTransf+1) ){
	  OutputTransformFileName  = argv[argcOutputTransf];
	  OutputTransform.open( OutputTransformFileName,
													std::ofstream::out | std::ofstream::trunc );
  }
  if(argc >= (argcOutputInvTransf+1) ){
	  OutputInvTransformFileName = argv[argcOutputInvTransf];
	  OutputInvTransform.open( OutputInvTransformFileName,
													std::ofstream::out | std::ofstream::trunc );
  }
	if( argc >= (argcInputTransfomX+1 ) ){
		inputTransformX = atof( argv[argcInputTransfomX] );
		inputTransformY = atof( argv[argcInputTransfomY] );
		inputTransformZ = atof( argv[argcInputTransfomZ] );
	}

  std::cout << std::endl << "Registration with moving "
		<< argv[argcMovingImage] << "\nand fixed\n" << argv[argcFixedImage]
		<< "\n"
		<< "Output File: " << argv[argcOutputImage] << std::endl;


////////////////////////////////////////////
//Definitions and instantiations
////////////////////////////////////////////

  const    unsigned int    ImageDimension = 3;
  typedef  short          PixelType;
  typedef itk::Image< unsigned short, 3 > ImageType;
  typedef itk::Image< unsigned short, 3 > OutputImageType;
  typedef itk::Image< signed short, 3> SignedImageType;
  typedef itk::Image<float, 3> DoubleImageType;
  typedef double CoordinateRepType;
  typedef itk::VersorRigid3DTransform< double > RTransformType;
  typedef itk::VersorRigid3DTransformOptimizer  ROptimizerType;
  typedef itk::LinearInterpolateImageFunction< SignedImageType, double>
																								InterpolatorType;
  typedef CommandIterationUpdate< ROptimizerType > RObserver;
  typedef itk::ResampleImageFilter< SignedImageType, SignedImageType> 
																								ResampleFilterType;
  typedef itk::ResampleImageFilter< SignedImageType, ImageType>
																								FinalResampleFilterType;
  typedef RTransformType::VersorType						VersorType;
  typedef VersorType::VectorType								VectorType;
  typedef ROptimizerType::ScalesType						OptimizerScalesType;
  typedef itk::Image<float, 3>									DoubleImageType;
  typedef itk::NormalizeImageFilter<SignedImageType,DoubleImageType>
																								NormalizeFilterType;
  typedef itk::ImageRegistrationMethod< DoubleImageType, DoubleImageType >
																								RRegistrationType;
  typedef itk::CenteredTransformInitializer< 
										RTransformType, 
										DoubleImageType, 
										DoubleImageType >						RTransformInitializerType;
  typedef itk::MattesMutualInformationImageToImageMetric<
										DoubleImageType, 
										DoubleImageType >						RMetricType;	
  typedef itk::LinearInterpolateImageFunction< DoubleImageType, double>
																								DoubleInterpolatorType;
  typedef itk::MetaImageIO											ImageIOType;
  typedef itk::ImageFileReader< SignedImageType >  ReaderType;
	typedef itk::CastImageFilter< SignedImageType, OutputImageType> CasterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  typedef itk::ImageRegionIteratorWithIndex< SignedImageType > IteratorType;
	typedef itk::IdentityTransform<double, ImageDimension>
																								IdentityTransformType;
		typedef itk::RecursiveGaussianImageFilter<
												SignedImageType,
												SignedImageType > GaussianFilterType;


///////////////////////////////////
// Reading Input Images
///////////////////////////////////

  ImageIOType::Pointer imageIO = ImageIOType::New(); 
  ReaderType::Pointer fixedImageReader = ReaderType::New();
  ReaderType::Pointer movingImageReader = ReaderType::New();
  
  fixedImageReader->SetImageIO( imageIO );
  fixedImageReader->SetFileName( argv[argcFixedImage] );
  SignedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  try{
      fixedImageReader->Update();
	}catch (itk::ExceptionObject &ex){
		std::cerr << "Error while reading the fixed image with " <<
			argv[argcFixedImage] << std::endl;
      std::cout << ex << std::endl;
      return EXIT_FAILURE;
	}
  
  movingImageReader->SetImageIO( imageIO );
  movingImageReader->SetFileName( argv[argcMovingImage] );
  try{
      movingImageReader->Update();
	}catch (itk::ExceptionObject &ex){
		std::cerr << "Error while reading the moving image with " <<
			argv[argcMovingImage] << std::endl;
    std::cout << ex << std::endl;
    return EXIT_FAILURE;
	}


	typedef itk::HistogramMatchingImageFilter<
				SignedImageType,
				SignedImageType > MatchingFilterType;
	MatchingFilterType::Pointer matcher = MatchingFilterType::New();
	matcher->SetInput( movingImageReader->GetOutput() );
	matcher->SetReferenceImage( fixedImageReader->GetOutput() );
	matcher->SetNumberOfHistogramLevels( 2048 );
	matcher->SetNumberOfMatchPoints( 5000 );
	//matcher->ThresholdAtMeanIntensityOn();
	matcher->Update();
  SignedImageType::Pointer movingImage = matcher->GetOutput();
	movingImage->DisconnectPipeline();

////////////////////////////////////

  ImageType::RegionType fixedRegion    =
																	fixedImage->GetLargestPossibleRegion();
  ImageType::SizeType fixedImageSize   = fixedRegion.GetSize();
  ImageType::IndexType fixedImageStart = fixedRegion.GetIndex();
  ImageType::RegionType movingRegion   = 
																	movingImage->GetLargestPossibleRegion();
  ImageType::SizeType movingImageSize  = movingRegion.GetSize();
  ImageType::IndexType movingImageStart = movingRegion.GetIndex();
  ImageType::SpacingType fixedSpacing  = fixedImage->GetSpacing();
  ImageType::SpacingType movingSpacing = movingImage->GetSpacing();

	double minMovingSpacing = movingSpacing[0];
	for(unsigned int d = 1; d< ImageDimension; d++){
		minMovingSpacing = vnl_math_min(minMovingSpacing,movingSpacing[d]);
	}
	double minFixedSpacing = fixedSpacing[0];
	for(unsigned int d = 1; d< ImageDimension; d++){
		minFixedSpacing = vnl_math_min(minFixedSpacing,fixedSpacing[d]);
	}

	double minSpacing = vnl_math_min(minFixedSpacing,minMovingSpacing);

	double spacingRatioFixed[ImageDimension];
	for(unsigned int d = 0; d< ImageDimension; d++){
		spacingRatioFixed[d] = fixedSpacing[d]/minSpacing;
	}
	double spacingRatioMoving[ImageDimension];
	for(unsigned int d = 0; d< ImageDimension; d++){
		spacingRatioMoving[d] = movingSpacing[d]/minSpacing;
	}

  std::cout << "Fixed image size: " << fixedImageSize << "\n"
						<< "Moving image size: " << movingImageSize << "\n"
						<< "Fixed image origin: " << fixedImage->GetOrigin() << "\n"
						<< "Moving image origin: " << movingImage->GetOrigin() << "\n" 
						<< "Fixed image spacing: " << fixedSpacing << "\n" 
						<< "Moving image spacing: " << movingSpacing << std::endl;

// Rigid registration
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	unsigned int it;
	unsigned int minlevel = 0;
	unsigned int numberOfIterations = atoi(argv[argcNbrOfIterations]);
	unsigned int nlevels = atoi(argv[argcNbrOfLevels]);

	//minimum size must be 2 in each dimension for Fixed;
	while(
		((pow(2,static_cast<float>(nlevels)))/spacingRatioFixed[0]>fixedImageSize[0]) |
		((pow(2,static_cast<float>(nlevels)))/spacingRatioFixed[1]>fixedImageSize[1]) |
		((pow(2,static_cast<float>(nlevels)))/spacingRatioFixed[2]>fixedImageSize[2])
	){
		nlevels=nlevels-1;
	}
	//minimum size must be 2 in each dimension for Moving;
	while(
		((pow(2,static_cast<float>(nlevels)))/spacingRatioMoving[0]>movingImageSize[0]) | 
		((pow(2,static_cast<float>(nlevels)))/spacingRatioMoving[1]>movingImageSize[1]) |
		((pow(2,static_cast<float>(nlevels)))/spacingRatioMoving[2]>movingImageSize[2])
	){
		nlevels=nlevels-1;
	}

	std::cout << "Number of Level: " << nlevels << std::endl;

	float maxNumberVoxels = 4e7;

	while(fixedImageSize[0]*fixedImageSize[1]*fixedImageSize[2]*spacingRatioFixed[0]*spacingRatioFixed[1]*spacingRatioFixed[2]/pow(2,3*static_cast<float>(minlevel)) > maxNumberVoxels) {
		minlevel=minlevel+1;
	}
	//maximum total size must be under 1Mvoxel
	while(movingImageSize[0]*movingImageSize[1]*movingImageSize[2]*spacingRatioMoving[0]*spacingRatioMoving[1]*spacingRatioMoving[2]/pow(2,3*static_cast<float>(minlevel)) > maxNumberVoxels) {
		minlevel=minlevel+1;
	}	

	if(minlevel>=nlevels){
		minlevel=nlevels-1;
	}


RTransformType::Pointer  Rtransform = RTransformType::New();

itk::Index<SignedImageType::ImageDimension> fixedCenterIndex =
									{{vnl_math_rnd(fixedImageSize[0]/2.0),
									vnl_math_rnd(fixedImageSize[1]/2.0),
									vnl_math_rnd(fixedImageSize[2]/2.0)}};
	itk::Index<SignedImageType::ImageDimension> movingCenterIndex =
									{{vnl_math_rnd(movingImageSize[0]/2.0),
									vnl_math_rnd(movingImageSize[1]/2.0),
									vnl_math_rnd(movingImageSize[2]/2.0)}};
	SignedImageType::PointType centerFixed;
 	fixedImage->TransformIndexToPhysicalPoint(fixedCenterIndex,centerFixed);
	SignedImageType::PointType centerMoving;
 	movingImage->TransformIndexToPhysicalPoint(movingCenterIndex,centerMoving);

	Rtransform->SetCenter( centerFixed );
	// set the ninitial translation: each center should coincide
	typedef itk::Vector<double, 3> VectorType;
  VectorType initTranslation;
	VectorType inputInitTranslation;
	inputInitTranslation[0] = inputTransformX;
	inputInitTranslation[1] = inputTransformY;
	inputInitTranslation[2] = inputTransformZ;

	for( unsigned int d=0; d<ImageDimension; d++ ){
		initTranslation[d] = centerMoving[d] - centerFixed[d] 
														+ inputInitTranslation[d] ;
	}
	Rtransform->SetTranslation( initTranslation );// centerMoving - centerFixed );

	VersorType     rotation;
	VectorType     axis;

	axis[0] = 0.0;
	axis[1] = 0.0;
	axis[2] = 1.0;
	const double angle = 0;

	rotation.Set(  axis, angle  );
	Rtransform->SetRotation( rotation );
	InterpolatorType::Pointer interpolator  = InterpolatorType::New();
	DoubleInterpolatorType::Pointer doubleInterpolator = 
																						DoubleInterpolatorType::New();
	ROptimizerType::ParametersType RfinalParameters =
																							Rtransform->GetParameters();

	std::cout << "Initial translation provided by user: ["
		<< inputTransformX << ", " << inputTransformY << ", "
		<< inputTransformZ <<"]" << std::endl;
	std::cout << "Initial translation: " << initTranslation
		<< " Init rotation: " << angle << "\n" << std::endl;

	RTransformInitializerType::Pointer Rinitializer =
		RTransformInitializerType::New();
	Rinitializer->SetTransform( Rtransform );
	//std::cout << "Transformation initial: "<< Rinitializer <<std::endl;

	////////////////////////////
	// For each level
	////////////////////////////
	for (it=nlevels;it>minlevel;it--){

		float movingFactors[ImageDimension];
		float fixedFactors[ImageDimension];

		for (int d = 0; d<ImageDimension; d++){
			fixedFactors[d]=(pow(2,static_cast<float>(it-1)))/spacingRatioFixed[d];
			if (fixedFactors[d]<1){
				fixedFactors[d]=1;
			}
			movingFactors[d]=(pow(2,static_cast<float>(it-1)))/spacingRatioMoving[d];
			if (movingFactors[d]<1){
				movingFactors[d]=1;
			}
		}

		// ***************** Smoothing and shrinking ************

		SignedImageType::Pointer tempFixedImage;

		// Fixed image
		GaussianFilterType::Pointer fixedSmootherX = GaussianFilterType::New();
		fixedSmootherX->SetInput( fixedImage );
		const double fixedSigmaX = fixedSpacing[0] * 0.5*( fixedFactors[0]);
		fixedSmootherX->SetSigma( fixedSigmaX );
		fixedSmootherX->SetDirection( 0 );
		fixedSmootherX->SetNormalizeAcrossScale( false );

		fixedSmootherX->Update();
		tempFixedImage = fixedSmootherX->GetOutput();
		tempFixedImage->DisconnectPipeline();
		
		GaussianFilterType::Pointer fixedSmootherY = GaussianFilterType::New();
		fixedSmootherY->SetInput( tempFixedImage );
		const double fixedSigmaY =fixedSpacing[1] * 0.5*( fixedFactors[1]);
		fixedSmootherY->SetSigma( fixedSigmaY );
		fixedSmootherY->SetDirection( 1 );
		fixedSmootherY->SetNormalizeAcrossScale( false );

		fixedSmootherY->Update();
		tempFixedImage->ReleaseData();
		tempFixedImage = fixedSmootherY->GetOutput();
		tempFixedImage->DisconnectPipeline();
		
		GaussianFilterType::Pointer fixedSmootherZ = GaussianFilterType::New();
		fixedSmootherZ->SetInput( tempFixedImage );
		const double fixedSigmaZ = fixedSpacing[2] * 0.5*(fixedFactors[2]);
		fixedSmootherZ->SetSigma( fixedSigmaZ );
		fixedSmootherZ->SetDirection( 2 );
		fixedSmootherZ->SetNormalizeAcrossScale( false );

		fixedSmootherZ->Update();
		tempFixedImage->ReleaseData();
		tempFixedImage = fixedSmootherZ->GetOutput();
		tempFixedImage->DisconnectPipeline();

		IdentityTransformType::Pointer identityTransform =IdentityTransformType::New();
		identityTransform->SetIdentity();

		ResampleFilterType::Pointer Fresampler = ResampleFilterType::New();
		Fresampler->SetTransform( identityTransform );
		Fresampler->SetInterpolator( interpolator);

		Fresampler->SetInput( tempFixedImage );
		typedef SignedImageType::SizeType::SizeValueType SizeValueType;
		SignedImageType::SizeType fSize;

		for (int d = 0; d<ImageDimension; d++){
			fSize[d] = static_cast<SizeValueType>(fixedImageSize[d]/fixedFactors[d]);
		}
		Fresampler->SetSize(fSize );
		Fresampler->SetOutputOrigin(  fixedImage->GetOrigin() );
		SignedImageType::SpacingType defSpacing;

		for (int d = 0; d<ImageDimension; d++){
			defSpacing[d] = fixedSpacing[d]*fixedFactors[d];
		}
		Fresampler->SetOutputSpacing( defSpacing );
		Fresampler->SetDefaultPixelValue( 0 );
		Fresampler->SetOutputDirection(fixedImage->GetDirection() );

		Fresampler->Update();
		tempFixedImage->ReleaseData();
		tempFixedImage = Fresampler->GetOutput();
		tempFixedImage->DisconnectPipeline();

		NormalizeFilterType::Pointer Fnormalizer = NormalizeFilterType::New();
		Fnormalizer->SetInput(tempFixedImage);
			
		Fnormalizer->Update();
		tempFixedImage->ReleaseData();
		DoubleImageType::Pointer doubleTempFixedImage = Fnormalizer->GetOutput();
		doubleTempFixedImage->DisconnectPipeline();
		
		// moving image
		SignedImageType::Pointer tempMovingImage;

		GaussianFilterType::Pointer movingSmootherX = GaussianFilterType::New();
		movingSmootherX->SetInput( movingImage );
		const double movingSigmaX = movingSpacing[0] * 0.5*( movingFactors[0]);
		movingSmootherX->SetSigma( movingSigmaX );
		movingSmootherX->SetDirection( 0 );
		movingSmootherX->SetNormalizeAcrossScale( false );

		movingSmootherX->Update();
		tempMovingImage = movingSmootherX->GetOutput();
		tempMovingImage->DisconnectPipeline();
		
		GaussianFilterType::Pointer movingSmootherY = GaussianFilterType::New();
		movingSmootherY->SetInput( tempMovingImage );
		const double movingSigmaY = movingSpacing[1] * 0.5*( movingFactors[1]);
		movingSmootherY->SetSigma( movingSigmaY );
		movingSmootherY->SetDirection( 1 );
		movingSmootherY->SetNormalizeAcrossScale( false );

		movingSmootherY->Update();
		tempMovingImage->ReleaseData();
		tempMovingImage = movingSmootherY->GetOutput();
		tempMovingImage->DisconnectPipeline();
		
		GaussianFilterType::Pointer movingSmootherZ = GaussianFilterType::New();
		movingSmootherZ->SetInput( tempMovingImage );
		const double movingSigmaZ = movingSpacing[2] * 0.5*( movingFactors[2]);
		movingSmootherZ->SetSigma( movingSigmaZ );
		movingSmootherZ->SetDirection( 2 );
		movingSmootherZ->SetNormalizeAcrossScale( false );

		movingSmootherZ->Update();
		tempMovingImage->ReleaseData();
		tempMovingImage = movingSmootherZ->GetOutput();
		tempMovingImage->DisconnectPipeline();
	
		ResampleFilterType::Pointer Mresampler = ResampleFilterType::New();
		Mresampler->SetTransform( identityTransform );
		Mresampler->SetInterpolator( interpolator);
		Mresampler->SetInput( tempMovingImage );				

		SignedImageType::SizeType mSize;

		for (int d = 0; d<ImageDimension; d++){
			mSize[d] = static_cast< SizeValueType>(movingImageSize[d]/movingFactors[d]);
		}
		std::cout << "sizes on the level:\nmoving: "
			<< mSize <<"\nfixed: " << fSize << std::endl;
		Mresampler->SetSize(mSize );
		Mresampler->SetOutputOrigin(  movingImage->GetOrigin() );
		Mresampler->SetOutputDirection(movingImage->GetDirection() );
		SignedImageType::SpacingType defSpacing2;

		for (int d = 0; d<ImageDimension; d++){
			defSpacing2[d] = movingSpacing[d]*movingFactors[d];
		}
		Mresampler->SetOutputSpacing( defSpacing2 );
		Mresampler->SetDefaultPixelValue( 0 );

		Mresampler->Update();
		tempMovingImage->ReleaseData();
		tempMovingImage = Mresampler->GetOutput();
		tempMovingImage->DisconnectPipeline();

		NormalizeFilterType::Pointer Mnormalizer = NormalizeFilterType::New();
		Mnormalizer->SetInput(tempMovingImage);
		
		Mnormalizer->Update();
		tempMovingImage->ReleaseData();
		DoubleImageType::Pointer doubleTempMovingImage = Mnormalizer->GetOutput();
		doubleTempMovingImage->DisconnectPipeline();


		// **************** Registration ***********

		Rtransform->SetParameters( RfinalParameters );

		RMetricType::Pointer  Rmetric = RMetricType::New();
		Rmetric -> ReinitializeSeed(76926294);
		double sample_ratio = 0.50;


		unsigned long numberOfSamples=fixedImageSize[0]*fixedImageSize[1]*fixedImageSize[2]*sample_ratio;
			Rmetric -> SetNumberOfSpatialSamples(numberOfSamples);

		ROptimizerType::Pointer Roptimizer = ROptimizerType::New();
		RRegistrationType::Pointer Rregistration = RRegistrationType::New();

		Rregistration->SetMetric( Rmetric );
		Rregistration->SetOptimizer( Roptimizer );
		Rregistration->SetInterpolator( doubleInterpolator );

		Rregistration->SetTransform( Rtransform );
		Rregistration->SetFixedImage( doubleTempFixedImage );
		Rregistration->SetMovingImage( doubleTempMovingImage );
		Rregistration->SetFixedImageRegion( fixedRegion );

		Rtransform->SetParameters( RfinalParameters );
		Rregistration->SetInitialTransformParameters( Rtransform->GetParameters() );

		OptimizerScalesType optimizerScales( Rtransform->GetNumberOfParameters() );

		//Scale so that the value will be close to 1
		//Versor will have a typical maximum value at 10 degrees
		//necessary if one of the directions does not have the same scale 
		optimizerScales[0] = 1.0*fixedSigmaX; 
		//necessary if one of the directions does not have the same scale 
		optimizerScales[1] = 1.0*fixedSigmaY; //made a change, it was fixedSigmaX
		optimizerScales[2] = 1.0*fixedSigmaZ; //made a change, it was fixedSigmaX

		//TODO: change that for application to GP
		//Translation will typically be about 50 mm x scale 
		//(if scale of 4 this gives about 20 cm)
		// [software guide] The optimizer scales the metrics (the gradient 
		// in this case) by the scales during each iteration.
		// Hence a large value of the center scale will prevent movement along 
		// the center during optimization.
		// WARNING: We want to constraint the translation along Z axis
		const double translationScale = 1 / (1e4) ; // 1 / (2e4);
		optimizerScales[3] = translationScale;
		optimizerScales[4] = translationScale;
		optimizerScales[5] = translationScale*100;
		Roptimizer->SetScales( optimizerScales );

		//maximum step length is about 1/20 of the maximum value set for
		//the parameters		
		Roptimizer->SetMaximumStepLength( 0.5 );
		Roptimizer->SetMinimumStepLength(0.0001);
		Roptimizer->SetNumberOfIterations( numberOfIterations*it );
		RObserver::Pointer Robserver = RObserver::New();
		Roptimizer->AddObserver( itk::IterationEvent(), Robserver );
		 		 
		try {
			Rregistration->Update();//StartRegistration();
		}catch( itk::ExceptionObject & err ){
			std::cerr << "ExceptionObject caught during registration!"
				<< std::endl;
			std::cerr << err << std::endl;
			exit( EXIT_FAILURE );
			//throw(-1);
		}

		RfinalParameters = Rregistration->GetLastTransformParameters();
		const unsigned int numberOfIterations=Roptimizer->GetCurrentIteration();
		const double bestValue = Roptimizer->GetValue();

		std::cout << std::endl ;
		std::cout << "Results : " << std::endl;
		std::cout << "Iterations    = " << numberOfIterations << std::endl;
		std::cout << "Metric value  = " << bestValue          << std::endl;

		Rtransform->SetParameters( RfinalParameters );
		RTransformType::MatrixType matrix = Rtransform->GetMatrix();
		RTransformType::OffsetType offset = Rtransform->GetOffset();
		RTransformType::CenterType center = Rtransform->GetCenter();

		std::cout << "Matrix = " << std::endl << matrix;
		std::cout << "Offset = " << std::endl << offset << std::endl;
		std::cout << "Center of rotation = " << std::endl << center 
		<< "\n" << std::endl;
	}//end for each iteration



  RTransformType::Pointer finalTransform = RTransformType::New();
  finalTransform->SetCenter( Rtransform->GetCenter() );
  finalTransform->SetParameters( RfinalParameters );

  ResampleFilterType::Pointer Rresampler = ResampleFilterType::New();
  Rresampler->SetTransform( finalTransform );
  Rresampler->SetInterpolator( interpolator);
  Rresampler->SetInput( movingImageReader->GetOutput() );
  Rresampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  Rresampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  Rresampler->SetOutputSpacing( fixedSpacing );
  Rresampler->SetDefaultPixelValue( 0 );

  try{
    Rresampler->Update();
	}catch ( itk::ExceptionObject & err ){
    std::cerr << "An error occured while resampling the moving image after rigid registration" << std::endl;
    return EXIT_FAILURE;
	}


// Write the output image
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  itk::TransformFileWriterTemplate<double>::Pointer inverseTransformWriter =
    itk::TransformFileWriterTemplate<double>::New();
  itk::TransformFileWriterTemplate<double>::Pointer transformWriter =
    itk::TransformFileWriterTemplate<double>::New();
#else
  itk::TransformFileWriter::Pointer inverseTransformWriter =
		itk::TransformFileWriter::New();
  itk::TransformFileWriter::Pointer transformWriter =
		itk::TransformFileWriter::New();
#endif

  if(argc>=7){  
	  std::cout << "Exporting transformation parameters to "
			<< argv[argcOutputTransf] << std::endl;
		transformWriter->SetFileName( argv[argcOutputTransf] );
		transformWriter->SetInput( finalTransform );
		try{
				transformWriter->Update();
		}catch(itk::ExceptionObject &err ){
			std::cerr << "Error while writing the inverse transform as: ";
			std::cerr << argv[argcOutputInvTransf] <<std::endl;
			std::cerr << err << std::endl;
		}
		/*
	  RTransformType::MatrixType matrix = finalTransform->GetMatrix();
	  RTransformType::OffsetType offset = finalTransform->GetOffset();
      OutputTransform << "Translation" << std::endl << offset << std::endl << std::endl;
      OutputTransform << "Rotation" << std::endl << matrix;
		*/
  }
  if(argc==8){  
	  std::cout << "Exporting inverse transformation parameters to "
			<< argv[argcOutputInvTransf] << std::endl;
		inverseTransformWriter->SetFileName( argv[argcOutputInvTransf] );
		inverseTransformWriter->SetInput(finalTransform->GetInverseTransform());
		try{
				inverseTransformWriter->Update();
		}catch(itk::ExceptionObject &err ){
			std::cerr << "Error while writing the inverse transform as: ";
			std::cerr << argv[argcOutputInvTransf] <<std::endl;
			std::cerr << err << std::endl;
			return(EXIT_FAILURE);
		}
	  /*
		 RTransformType::MatrixType matrix = finalTransform->GetInverseTransform()->GetMatrix();
	  RTransformType::OffsetType offset = finalTransform->GetInverseTransform()->GetOffset();
      OutputInvTransform << "Translation" << std::endl << offset << std::endl << std::endl;
      OutputInvTransform << "Rotation" << std::endl << matrix;
			*/
  }

   IteratorType itr( fixedImage, fixedImage->GetLargestPossibleRegion() );
   itr.GoToBegin();
   while ( ! itr.IsAtEnd() ){
     const ImageType::IndexType index = itr.GetIndex();
     PixelType v = (Rresampler->GetOutput())->GetPixel(index);
     itr.Set( v );
     ++itr;
   }

  std::cout << "Writing output image in "<< argv[argcOutputImage]
		<< std::endl;    
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput( fixedImage );
    WriterType::Pointer outputWriter = WriterType::New();
    std::string outputname = argv[argcOutputImage];
    outputWriter->SetFileName( outputname.c_str() );
    outputWriter->SetInput( caster->GetOutput() );
	try
    {
    outputWriter->Update();
	}
  catch( itk::ExceptionObject & excp )
    {
    std::cerr << "Exception thrown while writing the output image " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
    }


  std::cout << std::endl << "END   OF   THE   REGISTRATION   FROM   " << argv[2] << "   TO   " << argv[1] << std::endl << std::endl << std::endl;

  return 0;
}
