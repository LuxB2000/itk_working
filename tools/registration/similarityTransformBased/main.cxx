/*
 * A rigid registraiton method that takes into consideration the scaling
 */


// my libraries
#include "../../dcmreader/dcmreader.h"
#include "../../dcmwriter/dcmwriter.h"
#include "../../console_tools/color.h"

#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkTransformFileWriter.h"
#include "itkHistogramMatchingImageFilter.h"

// transform
#include "itkSimilarity3DTransform.h"
#include "itkCenteredTransformInitializer.h"
#include <itkIdentityTransform.h>

// optimizer
#include "itkVersorTransformOptimizer.h"
#include "itkImageRegistrationMethod.h"

// Metric
#include "itkMattesMutualInformationImageToImageMetric.h"

// smoothing
#include "itkDiscreteGaussianImageFilter.h"
#include "itkNormalizeImageFilter.h"
#include "itkResampleImageFilter.h"

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

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 4 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " FixedImagePath MovingImagePath outputImagePath"
			" [NbrOfLevels, default=2] [MinimalNbrOfIteration, default=10]"
			" [SavedDirectTransform, default=None]"
			" [SavedInverseTransform, default=None]"
			" [XInputTransform YInputTransform ZInputTransform]"<<std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argcFixedImage = 1, argcMovingImage = 2, argcOutputImage = 3,
							 argcNbrOfLevels = 4, argcNbrOfIterations = 5, 
							 argcDirectTransf = 6, argcInverTrans = 7, 
							 argcInputTransformX = 8, argcInputTransformY = 9, 
							 argcInputTransformZ = 10;

	std::string fixedImagePath = argv[argcFixedImage],
							movingImagePath = argv[argcMovingImage],
							outputPath = argv[argcOutputImage];

	// Default value
	unsigned int nbrOfIterations = 10;// atoi(argv[argcNbrOfIterations]);
	unsigned int nlevels = 2; // atoi(argv[argcNbrOfLevels]);
	float inputTransformX = 0, inputTransformY = 0, inputTransformZ = 0;
	std::string dirTransformPath = "", invTransformPath = "";

	if( argc > argcNbrOfLevels ){
		nlevels = atoi( argv[argcNbrOfLevels] );
	}
	if( argc > argcNbrOfIterations ){
		nbrOfIterations = atoi( argv[argcNbrOfIterations] );
	}
	if( argc > argcDirectTransf ){
		dirTransformPath = argv[argcDirectTransf];
	}
	if( argc > argcInverTrans ){
		invTransformPath = argv[argcInverTrans];
	}
	if( argc > argcInputTransformZ ){
		inputTransformX = atof( argv[argcInputTransformX] );
		inputTransformY = atof( argv[argcInputTransformY] );
		inputTransformZ = atof( argv[argcInputTransformZ] );
	}

	// plot inputs
	std::string msg;
	msg = std::string("Input parameters: ")
		+ std::string("\n- fixed image: ") + fixedImagePath
		+ std::string("\n- moving image: ") + movingImagePath
		+ std::string("\n- output image: ") + outputPath
		+ std::string("\n- ") + SSTR(nlevels) +  std::string(" levels")
		+ std::string("\n- ") + 
									SSTR(nbrOfIterations) +  std::string(" iterations")
		+ std::string("\n- ") +
									std::string("save direct transform ") + dirTransformPath
		+ std::string("\n- ") +
									std::string("save inverse transform ") + invTransformPath
		+ std::string("\n- ") +  std::string("inputTransform [")
						+ SSTR(inputTransformX) + std::string(", ") 
						+ SSTR(inputTransformY) +std::string(", ") +SSTR(inputTransformZ)
						+ std::string("]\n");
	blueMessage(msg);


	//
	// main type def
	//
  const    unsigned int    Dimension = 3;
  typedef  short          PixelType;
	// image types
  typedef itk::Image< PixelType, Dimension >							ImageType;
  typedef itk::Image< float, Dimension>										DoubleImageType;
	typedef DcmReader<ImageType, Dimension>									MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>									MyDCMWriterType;
  typedef itk::Similarity3DTransform< double >						RTransformType;
  typedef RTransformType::VersorType											VersorType;
	typedef itk::VersorTransformOptimizer										ROptimizerType;
  //typedef itk::VersorRigid3DTransformOptimizer  ROptimizerType;
  typedef itk::LinearInterpolateImageFunction< ImageType, double>
																													InterpolatorType;
  typedef itk::LinearInterpolateImageFunction< DoubleImageType, double>
																											DoubleInterpolatorType;
  typedef itk::CenteredTransformInitializer< 
										RTransformType, 
										DoubleImageType, 
										DoubleImageType >							RTransformInitializerType;
	typedef itk::RecursiveGaussianImageFilter<
												ImageType,
												ImageType >											GaussianFilterType;
  typedef itk::NormalizeImageFilter<ImageType,DoubleImageType>
																								NormalizeFilterType;
	typedef itk::IdentityTransform<double, Dimension>
																											IdentityTransformType;
  typedef itk::ResampleImageFilter< ImageType, ImageType> 
																													ResampleFilterType;
  typedef itk::MattesMutualInformationImageToImageMetric<
										DoubleImageType, 
										DoubleImageType >														RMetricType;
  typedef itk::ImageRegistrationMethod< DoubleImageType, DoubleImageType >
																								RRegistrationType;
  typedef ROptimizerType::ScalesType						OptimizerScalesType;
  typedef CommandIterationUpdate< ROptimizerType > RObserver;

	ImageType::Pointer fixedImage, movingImage;

	// read the images
	// and pre-processing
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( fixedImagePath.c_str() );
	try{
		dcmReader.Update();
		fixedImage = dcmReader.GetOutput();// dcmReader->GetOutput();
		fixedImage->DisconnectPipeline();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the fixed image with "
			<< fixedImagePath << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}
	MyDCMReaderType dcmReaderMoving = MyDCMReaderType();
	dcmReaderMoving.SetFileName( movingImagePath.c_str() );
	try{
		dcmReaderMoving.Update();
		movingImage = dcmReaderMoving.GetOutput();// dcmReader->GetOutput();
		movingImage->DisconnectPipeline();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the moving image with "
			<< fixedImagePath << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}
	typedef itk::HistogramMatchingImageFilter<
				ImageType,
				ImageType > MatchingFilterType;
	MatchingFilterType::Pointer matcher = MatchingFilterType::New();
	matcher->SetInput( movingImage );
	matcher->SetReferenceImage( fixedImage );
	matcher->SetNumberOfHistogramLevels( 2048 );
	matcher->SetNumberOfMatchPoints( 5000 );
	//matcher->ThresholdAtMeanIntensityOn();
	matcher->Update();
  movingImage = matcher->GetOutput();
	movingImage->DisconnectPipeline();

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
	for(unsigned int d = 1; d< Dimension; d++){
		minMovingSpacing = vnl_math_min(minMovingSpacing,movingSpacing[d]);
	}
	double minFixedSpacing = fixedSpacing[0];
	for(unsigned int d = 1; d< Dimension; d++){
		minFixedSpacing = vnl_math_min(minFixedSpacing,fixedSpacing[d]);
	}

	double minSpacing = vnl_math_min(minFixedSpacing,minMovingSpacing);

	double spacingRatioFixed[Dimension];
	for(unsigned int d = 0; d< Dimension; d++){
		spacingRatioFixed[d] = fixedSpacing[d]/minSpacing;
	}
	double spacingRatioMoving[Dimension];
	for(unsigned int d = 0; d< Dimension; d++){
		spacingRatioMoving[d] = movingSpacing[d]/minSpacing;
	}

  std::cout << "Fixed image size: " << fixedImageSize << "\n"
						<< "Moving image size: " << movingImageSize << "\n"
						<< "Fixed image origin: " << fixedImage->GetOrigin() << "\n"
						<< "Moving image origin: " << movingImage->GetOrigin() << "\n" 
						<< "Fixed image spacing: " << fixedSpacing << "\n" 
						<< "Moving image spacing: " << movingSpacing << std::endl;

	//
	// Perform registration
	//
	unsigned int it=0;
	unsigned int minlevel = 0;
	
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

	itk::Index<ImageType::ImageDimension> fixedCenterIndex =
									{{vnl_math_rnd(fixedImageSize[0]/2.0),
									vnl_math_rnd(fixedImageSize[1]/2.0),
									vnl_math_rnd(fixedImageSize[2]/2.0)}};
	itk::Index<ImageType::ImageDimension> movingCenterIndex =
									{{vnl_math_rnd(movingImageSize[0]/2.0),
									vnl_math_rnd(movingImageSize[1]/2.0),
									vnl_math_rnd(movingImageSize[2]/2.0)}};
	ImageType::PointType centerFixed;
 	fixedImage->TransformIndexToPhysicalPoint(fixedCenterIndex,centerFixed);
	ImageType::PointType centerMoving;
 	movingImage->TransformIndexToPhysicalPoint(movingCenterIndex,centerMoving);
	Rtransform->SetCenter( centerFixed );

	// set the ninitial translation: each center should coincide
	typedef itk::Vector<double, 3> VectorType;
  VectorType initTranslation;
	VectorType inputInitTranslation;
	inputInitTranslation[0] = inputTransformX;
	inputInitTranslation[1] = inputTransformY;
	inputInitTranslation[2] = inputTransformZ;

	for( unsigned int d=0; d<Dimension; d++ ){
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

		float movingFactors[Dimension];
		float fixedFactors[Dimension];

		for (int d = 0; d<Dimension; d++){
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
		ImageType::Pointer tempFixedImage;

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
		typedef ImageType::SizeType::SizeValueType SizeValueType;
		ImageType::SizeType fSize;

		for (int d = 0; d<Dimension; d++){
			fSize[d] = static_cast<SizeValueType>(fixedImageSize[d]/fixedFactors[d]);
		}
		Fresampler->SetSize(fSize );
		Fresampler->SetOutputOrigin(  fixedImage->GetOrigin() );
		ImageType::SpacingType defSpacing;

		for (int d = 0; d<Dimension; d++){
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
		ImageType::Pointer tempMovingImage;

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

		ImageType::SizeType mSize;

		for (int d = 0; d<Dimension; d++){
			mSize[d] = static_cast< SizeValueType>(movingImageSize[d]/movingFactors[d]);
		}
		std::cout << "sizes on the level:\nmoving: "
			<< mSize <<"\nfixed: " << fSize << std::endl;
		Mresampler->SetSize(mSize );
		Mresampler->SetOutputOrigin(  movingImage->GetOrigin() );
		Mresampler->SetOutputDirection(movingImage->GetDirection() );
		ImageType::SpacingType defSpacing2;

		for (int d = 0; d<Dimension; d++){
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
		const double translationScale = 1 / (1e2) ; // 1 / (2e4);
		optimizerScales[3] = translationScale;
		optimizerScales[4] = translationScale;
		optimizerScales[5] = translationScale*100;
		// scale factor ... TODO: find a correct value
		optimizerScales[6] = 1;
		Roptimizer->SetScales( optimizerScales );

		//maximum step length is about 1/20 of the maximum value set for
		//the parameters		
		Roptimizer->SetMaximumStepLength( 0.5 );
		Roptimizer->SetMinimumStepLength(0.0001);
		Roptimizer->SetNumberOfIterations( nbrOfIterations*it );
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
		RTransformType::ScaleType scale =  Rtransform->GetScale();

		std::cout << "Matrix = " << std::endl << matrix;
		std::cout << "Offset = " << std::endl << offset << std::endl;
		std::cout << "Center of rotation = " << center 
		<< "\nScale: " << scale
		<< "\n" << std::endl;

	}// end for each level

	//
	// Apply the deformations to moving
	//


  RTransformType::Pointer finalTransform = RTransformType::New();
  finalTransform->SetCenter( Rtransform->GetCenter() );
  finalTransform->SetParameters( RfinalParameters );

  ResampleFilterType::Pointer Rresampler = ResampleFilterType::New();
  Rresampler->SetTransform( finalTransform );
  Rresampler->SetInterpolator( interpolator);
  Rresampler->SetInput( dcmReaderMoving.GetOutput() );
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

  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
   IteratorType itr( fixedImage, fixedImage->GetLargestPossibleRegion() );
   itr.GoToBegin();
   while ( ! itr.IsAtEnd() ){
     const ImageType::IndexType index = itr.GetIndex();
     PixelType v = (Rresampler->GetOutput())->GetPixel(index);
     itr.Set( v );
     ++itr;
   }

	typedef itk::CastImageFilter< ImageType, ImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();
	caster->SetInput( fixedImage );

	//////////////////////////////////////////////////////////////////////////////
	// write the transform
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
  if(dirTransformPath.length()>0){
	  msg = std::string("Exporting transformation parameters to ") + 
			dirTransformPath;
		finalMessage( msg );
		transformWriter->SetFileName( dirTransformPath );
		transformWriter->SetInput( finalTransform );
		try{
				transformWriter->Update();
		}catch(itk::ExceptionObject &err ){
			std::cerr << "Error while writing the direct transform as: ";
			std::cerr << dirTransformPath <<std::endl;
			std::cerr << err << std::endl;
		}
  }
  if(invTransformPath.length()>0){
	  msg = std::string("Exporting inverse transformation parameters to ")
			+ invTransformPath;
		finalMessage(msg);
		inverseTransformWriter->SetFileName( invTransformPath );
		inverseTransformWriter->SetInput(finalTransform->GetInverseTransform());
		try{
				inverseTransformWriter->Update();
		}catch(itk::ExceptionObject &err ){
			std::cerr << "Error while writing the inverse transform as: ";
			std::cerr << invTransformPath <<std::endl;
			std::cerr << err << std::endl;
			return(EXIT_FAILURE);
		}
	}
	
	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( outputPath.c_str() );
	writer.SetInput( caster->GetOutput() );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the volume as: "
			<< outputPath << std::endl;
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}
	
	msg = "Output volume wrote with " + outputPath;
	finalMessage( msg );
	std::cout << std::endl;
}


