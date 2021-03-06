/*
 * Translation 3D Mutli-Resolution Normalized Registration. The registration is performed
 * between two volumes FIXED and MOVING and provided a 3D transformation M. 
 *
 * The registration is defined as:
 * transform			Translation
 * optimizer			GradientDescentOptimizer (for noisy images)
 * metric					Mutual Information (multi-modal images)
 * interpolator		LinearInterpolateImageFunction
 * registration		multiresolution registration
 *
 * INPUTS
 * ------
 *  * The fixed image path
 *  * The moving image path
 *  * The deformed image path 
 *	* The transform object path 
 *	* The transform inverse object path
 *	* Init translation in mm (optional, default 0 0 0 )
 *	* Ratio of pixel to take into account - between 1 and 100 (optional, default is 40)
 *	* Moment scale (optional, default 1)
 *	* Translation scale (optional, default 1/1e5)
 *
 * CONSOLE OUTPUTS
 * ---------------
 * (iter number)_[current metric]_[current value]
 * These are generated by CommandIterationUpdate
 * When the algo is finished, a final line is displayed informing what is the stop criterion.
 * Then the numerical results are display after \n\n
 *
 * OUTPUTS
 * -------
 * * The deformed volume
 * * The transform object and its inverse
 * Recommendation: the output volumes are written on the disks in /tmp/
 *
 * WARNINGS: The outputs are valid for UNIX station and they will be deleted at each reboot of the system.
 *
 *
 * EXAMPLES
 * --------
 * ./registration fixedImage.mha movingImage.mha /tmp/deformedImage.mha /tmp/transform-direct.tfm /tmp/transform-inverse.tfm
 * ./registration fixedImage.mha movingImage.mha /tmp/deformedImage.mha /tmp/transform.tfm /tmp/transform-inverse.tfm -15 0 -0.2 65 1 0.001
 *
 *
 *
 * Mostly inspirited by ImageRegistration2.cxx and
 * MultiResImageRegistration2.cxx from "ITK user guide"
 *
 * 2014
 * Author: J Plumat
 *
 */

// reader - writer
#include "../../dcmreader/dcmreader.h"
#include "../../dcmwriter/dcmwriter.h"

// affine transform object
#include "../lib/translationNormalisedMultiRegistration.h"


// general
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkVector.h"

// Normalizer and smoother
#include "itkNormalizeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

// measure time
#include <time.h>

// matrix computation
#include "itkVersor.h"

// viewer
#include "QuickView.h"

/*
 * MAIN
 */
int main( int argc, char *argv[] )
{
	std::cout<<" -- "<<argv[0]<<" -- "<<std::endl;
	// check parameters
  float momentScale = 1;
  float translationScale = 1.0 / 1; //1e5
	int ratio = 100;
	if( argc < 6 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImageFile  movingImageFile deformedImage";
    std::cerr << " transform inverseTransform\n"
			<< "Optional: [initTranslationX initTranslationY initTranslationZ] "
			<< "[pixelRatio=" << ratio << "]"
			<< "[momentScale=" << momentScale << "]"
			<< "[translationScale=" << translationScale << "]"<<
			std::endl;
    return EXIT_FAILURE;
    }

	const int argvFixedImageFile=1, argvMovingImageFile=2,
			argvOutputImageFile=3, argvOutputTransform=4,
			argvOutputInverseTransform=5, 
			argvTx=6, argvTy=7, argvTz=8, argvRatio=9, argvMomentScale=10,
			argvTranslationScale=11;

	const bool show = true;

	// starting time
	time_t beginTime,endTime;
	time(&beginTime);

	// default values
	itk::Vector<double, 3> initTranslation;
	initTranslation[0] = 0;
	initTranslation[1] = 0;
	initTranslation[2] = 0;
	if( argc >= 12 ) translationScale = atof(argv[argvTranslationScale]);
	if( argc >= 11 ) momentScale = atof(argv[argvMomentScale]);
	if( argc >= 10 ) ratio = atoi( argv[argvRatio] );
	if( argc >= 9 ){
		initTranslation[0] = (double) atof( argv[argvTx] );
		initTranslation[1] = (double) atof( argv[argvTy] );
		initTranslation[2] = (double) atof( argv[argvTz] );
	}
	
	std::cout << "Inputs: Tx,Ty,Tz= " << initTranslation
		<< " ratio= " << ratio << " translationScale= " << translationScale
		<< " momentScale= " << momentScale
		<< std::endl;


	//
	// CONST
	//
	const int nbrOfLevel = 3; // min size is along z and is 32 for GP
	const int maxIter = 200;
	//const double maxStepLength = 0.2000, minStepLength = 0.0001;

	//
	// main type def
	//
  const    unsigned int    Dimension = 3;
  typedef  float           PixelType;
  typedef  short					 OutputPixelType;
  typedef  short					 MaskPixelType;

	int defaultOutputPixel = 0;

  typedef itk::Image< PixelType, Dimension >        FixedImageType;
  typedef itk::Image< PixelType, Dimension >        MovingImageType;
  typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;

	MovingImageType::Pointer movingImage;
	FixedImageType::Pointer fixedImage, result;
	
	// transform object
	typedef AffineMultiRegistration< MovingImageType,
																		FixedImageType > TransformObjectType;

	//
	// Readers
	//
	typedef DcmReader<FixedImageType,Dimension> FixedDCMReaderType;
	FixedDCMReaderType* fixedDCMReader =
														new FixedDCMReaderType(argv[argvFixedImageFile]);
	typedef DcmReader<MovingImageType,Dimension> MovingDCMReaderType;
	FixedDCMReaderType* movingDCMReader =
													new MovingDCMReaderType(argv[argvMovingImageFile]);

	try{
		std::cout<<"Read the fixed image with path: ";
		std::cout<<argv[argvFixedImageFile]<<std::endl;
		fixedDCMReader->Update();
		fixedImage = fixedDCMReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr<<"Error while reading the fixed image."<<std::endl;
		std::cerr<<err<<std::endl;
		exit(EXIT_FAILURE);
	}

	try{
		std::cout<<"Read the moving image with path: ";
		std::cout<<argv[argvMovingImageFile]<<std::endl;
		movingDCMReader->Update();
		movingImage = movingDCMReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr<<"Error while reading the moving image."<<std::endl;
		std::cerr<<err<<std::endl;
		exit(EXIT_FAILURE);
	}

	// Normalizer and smoother
	typedef itk::NormalizeImageFilter< FixedImageType,FixedImageType>
																							FixedNormalizeFilterType;
	typedef itk::NormalizeImageFilter< MovingImageType,MovingImageType>
																							MovingNormalizeFilterImageType;

	typedef itk::DiscreteGaussianImageFilter< FixedImageType, FixedImageType> 
																							GaussianFilterType;

	FixedNormalizeFilterType::Pointer fixedNormalizer =
																				FixedNormalizeFilterType::New();
	MovingNormalizeFilterImageType::Pointer movingNormalizer =
																			MovingNormalizeFilterImageType::New();
	GaussianFilterType::Pointer fixedSmoother = GaussianFilterType::New();
	GaussianFilterType::Pointer movingSmoother = GaussianFilterType::New();

	fixedNormalizer->SetInput( fixedImage );
	movingNormalizer->SetInput( movingImage );

	fixedSmoother->SetVariance( 2.5 );
	fixedSmoother->SetInput( fixedNormalizer->GetOutput() );
	fixedSmoother->Update();
	movingSmoother->SetInput( movingNormalizer->GetOutput() );
	movingSmoother->SetVariance( 2.5 );
	movingSmoother->Update();

	// Registration

	typedef TransformObjectType::TransformType TransformType;

	TransformObjectType* doRegistration = new TransformObjectType();
	doRegistration->SetMovingImage( movingNormalizer->GetOutput() );
	doRegistration->SetFixedImage( fixedNormalizer->GetOutput() );
	doRegistration->SetMomentScale( momentScale );
	doRegistration->SetTranslationScale( translationScale );
	doRegistration->SetPixelRatio( ratio );
	doRegistration->SetNumberOfLevels( nbrOfLevel );
	doRegistration->SetInitTranslation( initTranslation[0],
			initTranslation[1],
			initTranslation[2]);
	try{
		doRegistration->Update();
	}catch(itk::ExceptionObject &err){
		std::cerr <<"Error during registration process: "<<std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	TransformObjectType::OptimizerType::ParametersType finalParameters =
			doRegistration->GetFinalParameters();
	TransformObjectType::OptimizerType::ParametersType fixedParameters = 
			doRegistration->GetFixedParameters();
	
	//
	// deformed the moving image and save the deformed.
	//
  typedef itk::ResampleImageFilter< 
                            MovingImageType, 
                            FixedImageType >       ResampleFilterType;
  typedef itk::CastImageFilter< 
                        MovingImageType,
                        OutputImageType >          CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;

  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  WriterType::Pointer         writer    =  WriterType::New();
  CastFilterType::Pointer     caster    =  CastFilterType::New();

	TransformType::Pointer finalTransform = TransformType::New();
  finalTransform->SetParameters( finalParameters );
/*			registration->GetLastTransformParameters() ); */

  finalTransform->SetFixedParameters( fixedParameters );
/*			transform->GetFixedParameters() ); */

	// set transform and initial moving image
  resampler->SetTransform( finalTransform );
  resampler->SetInput(     movingImage );

  resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  resampler->SetOutputSpacing( fixedImage->GetSpacing() );
  resampler->SetOutputDirection( fixedImage->GetDirection() );
  resampler->SetDefaultPixelValue( defaultOutputPixel );

	resampler->Update();
	result = resampler->GetOutput();

  writer->SetFileName( argv[argvOutputImageFile] );
  caster->SetInput( resampler->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
	
	// write the deformed moving image
	try
		{
			std::cout<<"Write the deformed image as: ";
			std::cout<<argv[argvOutputImageFile]<<std::endl;
			writer->Update();
		}
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught during the deformed";
		std::cerr << " moving image writing ! Proposed filename is ";
		std::cerr << argv[argvOutputImageFile]<<"." << std::endl; 
    std::cerr << err << std::endl; 
		return EXIT_FAILURE;
		}


	//
	// Write the transform and the inverse tranform
	//
	#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  itk::TransformFileWriterTemplate<double>::Pointer transformWriter =
    itk::TransformFileWriterTemplate<double>::New();
#else
  itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
#endif
	transformWriter->SetInput(finalTransform);
	transformWriter->SetFileName(argv[argvOutputTransform]);
	std::cout << "Write the direct transform as: ";
	std::cout << argv[argvOutputTransform] << std::endl;
	try
	{
		transformWriter->Update();
	}catch( itk::ExceptionObject &err)
	{
		std::cerr << "Error while writing the direct transform as: ";
		std::cerr << argv[argvOutputTransform] << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}


	#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  itk::TransformFileWriterTemplate<double>::Pointer inverseTransformWriter =
    itk::TransformFileWriterTemplate<double>::New();
#else
  itk::TransformFileWriter::Pointer inverseTransformWriter =
		itk::TransformFileWriter::New();
#endif
	inverseTransformWriter->SetInput(finalTransform->GetInverseTransform());
	inverseTransformWriter->SetFileName(argv[argvOutputInverseTransform]);
	std::cout << "Write the inverse transform as: ";
	std::cout << argv[argvOutputInverseTransform] << std::endl;
	try
	{
		inverseTransformWriter->Update();
	}catch(itk::ExceptionObject &err ){
		std::cerr << "Error while writing the inverse transform as: ";
		std::cerr << argv[argvOutputInverseTransform] <<std::endl;
		std::cerr << err << std::endl;
		return(EXIT_FAILURE);
	}

	/*
	if( show )
	{
		FixedImageType::IndexType start;
		start.Fill(0);
		FixedImageType::SizeType size;
		size.Fill
		QuickView viewer;
		viewer.AddImage( fixedImage.GetPointer(), true, "Fixed Image" );
		viewer.AddImage( movingImage.GetPointer(), true, "Moving Image" );
		viewer.AddImage( result.GetPointer(), true, "Result Image" );
	}
	*/

	std::cout << std::endl;
	exit(EXIT_SUCCESS);
}//end main
