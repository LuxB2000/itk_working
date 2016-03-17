#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

/*
 * Authors: Jonathan Orban de Xirvy and Guillaime Janssens, UCL
 *					Plumat Jerome, UoA
  std::ofstream OutputTransform;
 */

#include <stdlib.h>
#include <stdio.h>
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkShrinkImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkImage.h"
#include "itkAffineTransform.h"
#include "itkCenteredTransformInitializer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMetaImageIO.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkResampleImageFilter.h"
#include "itkCommand.h"
#include "itkSamplesMeanSquaresImageToImageMetric.h"
#include "itkSamplesNormalizedCorrelationImageToImageMetric.h"
//#include "itkMattesMutualInformationImageToImageMetricNew.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkTransformFileWriter.h"

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
	std::cout << " = " << (double)(optimizer->GetValue())<<std::endl;
    //std::cout << " = " << (1.0 + (double)(optimizer->GetValue()))*100.0<<std::endl;
}
};

// Main
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
  if( (argc<6) || (argc>9) ){
    std::cerr << "Wrong number of Parameters (" << argc << ")" << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " FixedImageDir MovingImageDir OutputImageDir Number_of_Levels Min_Number_of_Iterations (outputTransformFile)" << std::endl << std::endl;
    std::cerr << "Provided parameters :" << std::endl;
	for (int i = 1; i < argc; i++){
		std::cout << "ARG " << i << " : " << argv[i] << std::endl;}
    return 1;
    }

	const unsigned int argcDTransform = 6, argcITransform = 7;

  const char* OutputTransformFileName;
  std::cout << std::endl << "BEGIN  OF   THE   REGISTRATION   FROM   " << argv[2] << "   TO   " << argv[1] << std::endl;
  std::cout << "WRITING   RESULTS   IN   " << argv[3] << ".  (Number of params = " << argc << ")" << std::endl;

	bool saveDTransform = false, saveITransform = false;
  if(argc >= 7){	  
		std::cout << "Save Direct Transform to " << argv[argcDTransform] << std::endl;
		saveDTransform = true;
  }
	if( argc == 8 ){
		std::cout << "Save inverse transform to " << argv[argcITransform] << std::endl;
		saveITransform = true;
	}

//Definitions and instantiations
////////////////////////////////////////////

  const    unsigned int    ImageDimension = 3;
  typedef  unsigned short          PixelType;
  typedef itk::Image< PixelType, ImageDimension >  ImageType;
  typedef double CoordinateRepType;
  typedef itk::AffineTransform< double , ImageDimension> RTransformType;
  typedef itk::RegularStepGradientDescentOptimizer           ROptimizerType;
  typedef itk::SamplesMeanSquaresImageToImageMetric< ImageType, ImageType >  RMetricType;
  //typedef itk::SamplesNormalizedCorrelationImageToImageMetric< ImageType, ImageType >  RMetricType;
  //typedef itk::MattesMutualInformationImageToImageMetricNew< ImageType, ImageType >  RMetricType;
  typedef itk::LinearInterpolateImageFunction< ImageType, double>  InterpolatorType;
  typedef itk::ImageRegistrationMethod< ImageType, ImageType > RRegistrationType;
  typedef CommandIterationUpdate< ROptimizerType > RObserver;
  typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
  typedef itk::CenteredTransformInitializer< RTransformType, ImageType, ImageType   >  RTransformInitializerType;
  typedef ROptimizerType::ScalesType       OptimizerScalesType;
  typedef itk::MetaImageIO       ImageIOType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;


// Reading Input Images
///////////////////////////////////

  ImageIOType::Pointer imageIO = ImageIOType::New(); 
  ReaderType::Pointer fixedImageReader = ReaderType::New();
  ReaderType::Pointer movingImageReader = ReaderType::New();
  
  fixedImageReader->SetImageIO( imageIO );
  fixedImageReader->SetFileName( argv[1] );
  ImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  try{
      fixedImageReader->Update();}
    catch (itk::ExceptionObject &ex){
      std::cout << ex << std::endl;
      return EXIT_FAILURE;}   
  
  movingImageReader->SetImageIO( imageIO );
  movingImageReader->SetFileName( argv[2] );
  ImageType::Pointer movingImage = movingImageReader->GetOutput();
  try{
      movingImageReader->Update();}
    catch (itk::ExceptionObject &ex){
      std::cout << ex << std::endl;
      return EXIT_FAILURE;}


////////////////////////////////////

  ImageType::RegionType fixedRegion = fixedImage->GetLargestPossibleRegion();
  ImageType::SizeType fixedImageSize = fixedRegion.GetSize();
  ImageType::IndexType fixedImageStart = fixedRegion.GetIndex();
  ImageType::RegionType movingRegion = movingImage->GetLargestPossibleRegion();
  ImageType::SizeType movingImageSize = movingRegion.GetSize();
  ImageType::IndexType movingImageStart = movingRegion.GetIndex();
  ImageType::SpacingType fixedSpacing = fixedImage->GetSpacing();
  ImageType::SpacingType movingSpacing = movingImage->GetSpacing();

  std::cout << "Fixed image size: " << fixedImageSize << std::endl;
  std::cout << "Moving image size: " << movingImageSize << std::endl;
  std::cout << "Fixed image origin: " << fixedImage->GetOrigin() << std::endl;
  std::cout << "Moving image origin: " << movingImage->GetOrigin() << std::endl;
  std::cout << "Fixed image spacing: " << fixedSpacing << std::endl;
  std::cout << "Moving image spacing: " << movingSpacing << std::endl;

// Rigid registration
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int it;
int minlevel = 0;

int nlevel=atoi(argv[4]);

while(((vnl_math_sqr((pow(2,static_cast<float>(nlevel-1)))))>fixedImageSize[0]) | ((vnl_math_sqr((pow(2,static_cast<float>(nlevel-1)))))>fixedImageSize[1]) | ((vnl_math_sqr((pow(2,static_cast<float>(nlevel-1)))))>fixedImageSize[2])) {
nlevel=nlevel-1;
}
while(fixedImageSize[0]*fixedImageSize[1]*fixedImageSize[2]/pow(2,3*static_cast<float>(minlevel)) > 3e6) {
minlevel=minlevel+1;
}
if(minlevel>=nlevel){
minlevel=nlevel-1;
}
unsigned int factors[ImageDimension];

RTransformType::Pointer  Rtransform = RTransformType::New();
RTransformType::InputPointType centerFixed;
RTransformType::InputPointType centerMoving;

  centerFixed[0] = fixedImage->GetOrigin()[0] + fixedSpacing[0] * fixedImageSize[0] / 2.0;
  centerFixed[1] = fixedImage->GetOrigin()[1] + fixedSpacing[1] * fixedImageSize[1] / 2.0;
  centerFixed[2] = fixedImage->GetOrigin()[2] + fixedSpacing[2] * fixedImageSize[2] / 2.0;
  centerMoving[0] = movingImage->GetOrigin()[0] + movingSpacing[0] * movingImageSize[0] / 2.0;
  centerMoving[1] = movingImage->GetOrigin()[1] + movingSpacing[1] * movingImageSize[1] / 2.0;
  centerMoving[2] = movingImage->GetOrigin()[2] + movingSpacing[2] * movingImageSize[2] / 2.0;

Rtransform->SetCenter( centerFixed );
Rtransform->SetTranslation( centerMoving - centerFixed );

InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
ROptimizerType::ParametersType RfinalParameters = Rtransform->GetParameters();

RTransformInitializerType::Pointer Rinitializer =
                                          RTransformInitializerType::New();
Rinitializer->SetTransform( Rtransform );

for (it=nlevel;it>minlevel;it--)
{
  std::cout << std::endl << "--> Starting rigid registration level " << nlevel-it+1 << " on " << nlevel-minlevel << std::endl;

  factors[0]=static_cast<int>(pow(2,static_cast<float>(it-1)));
  factors[1]=static_cast<int>(pow(2,static_cast<float>(it-1)));
  factors[2]=static_cast<int>(pow(2,static_cast<float>(it-1)));

  // ***************** Smoothing and shrinking ************

  typedef itk::RecursiveGaussianImageFilter<ImageType,ImageType > GaussianFilterType;

  // Fixed image
  GaussianFilterType::Pointer fixedSmootherX = GaussianFilterType::New();
  GaussianFilterType::Pointer fixedSmootherY = GaussianFilterType::New();
  GaussianFilterType::Pointer fixedSmootherZ = GaussianFilterType::New();
  fixedSmootherX->SetInput( fixedImage );
  fixedSmootherY->SetInput( fixedSmootherX->GetOutput() );
  fixedSmootherZ->SetInput( fixedSmootherY->GetOutput() );
  const double fixedSigmaX = fixedSpacing[0] * vnl_math_sqr(0.5*static_cast<float>( factors[0]));
  const double fixedSigmaY = fixedSpacing[1] * vnl_math_sqr(0.5*static_cast<float>( factors[1]));
  const double fixedSigmaZ = fixedSpacing[2] * vnl_math_sqr(0.5*static_cast<float>( factors[2]));
  std::cout << "Setting variance for fixed image smoothing: " << fixedSigmaX<< " " << fixedSigmaY << " " << fixedSigmaZ << std::endl;
  fixedSmootherX->SetSigma( fixedSigmaX );
  fixedSmootherY->SetSigma( fixedSigmaY );
  fixedSmootherZ->SetSigma( fixedSigmaZ );
  fixedSmootherX->SetDirection( 0 );
  fixedSmootherY->SetDirection( 1 );
  fixedSmootherZ->SetDirection( 2 );
  fixedSmootherX->SetNormalizeAcrossScale( false );
  fixedSmootherY->SetNormalizeAcrossScale( false );
  fixedSmootherZ->SetNormalizeAcrossScale( false );

  typedef itk::ShrinkImageFilter<ImageType,ImageType> FixedShrinkType;
  FixedShrinkType::Pointer fixedShrink= FixedShrinkType::New();
  fixedShrink->SetInput(fixedSmootherZ->GetOutput());
  fixedShrink->SetShrinkFactors( factors );
  fixedShrink->Update();

  // Moving image
  if( (movingSpacing[0]>=2.0*fixedSpacing[0]) & (factors[0]>=2) ){
  factors[0]=static_cast<int>(factors[0]/2);}
  if( (movingSpacing[1]>=2.0*fixedSpacing[1]) & (factors[1]>=2) ){
  factors[1]=static_cast<int>(factors[1]/2);}
  if( (movingSpacing[2]>=2.0*fixedSpacing[2]) & (factors[2]>=2) ){
  factors[2]=static_cast<int>(factors[2]/2);}

  GaussianFilterType::Pointer movingSmootherX = GaussianFilterType::New();
  GaussianFilterType::Pointer movingSmootherY = GaussianFilterType::New();
  GaussianFilterType::Pointer movingSmootherZ = GaussianFilterType::New();
  movingSmootherX->SetInput( movingImage );
  movingSmootherY->SetInput( movingSmootherX->GetOutput() );
  movingSmootherZ->SetInput( movingSmootherY->GetOutput() );
  const double movingSigmaX = movingSpacing[0] * vnl_math_sqr(0.5*static_cast<float>( factors[0]));
  const double movingSigmaY = movingSpacing[1] * vnl_math_sqr(0.5*static_cast<float>( factors[1]));
  const double movingSigmaZ = movingSpacing[2] * vnl_math_sqr(0.5*static_cast<float>( factors[2]));
  std::cout << "Setting variance for moving image smoothing: " << movingSigmaX<< " " << movingSigmaY << " " << movingSigmaZ << std::endl;
  movingSmootherX->SetSigma( movingSigmaX );
  movingSmootherY->SetSigma( movingSigmaY );
  movingSmootherZ->SetSigma( movingSigmaZ );
  movingSmootherX->SetDirection( 0 );
  movingSmootherY->SetDirection( 1 );
  movingSmootherZ->SetDirection( 2 );
  movingSmootherX->SetNormalizeAcrossScale( false );
  movingSmootherY->SetNormalizeAcrossScale( false );
  movingSmootherZ->SetNormalizeAcrossScale( false );

  typedef itk::ShrinkImageFilter<ImageType,ImageType> MovingShrinkType;
  MovingShrinkType::Pointer movingShrink= MovingShrinkType::New();
  movingShrink->SetInput(movingSmootherZ->GetOutput());
  movingShrink->SetShrinkFactors( factors );
  movingShrink->Update();


  // **************** Registration ***********

  Rtransform->SetParameters( RfinalParameters );

  RMetricType::Pointer  Rmetric = RMetricType::New();
  Rmetric -> ReinitializeSeed(76926294);
  double sample_ratio = 0.15;

  Rmetric -> SetNumberOfSamples((unsigned long) ( sample_ratio*(double)((  ((fixedRegion.GetSize())[0]) * ((fixedRegion.GetSize())[1]) * ((fixedRegion.GetSize())[2]) ))));

  ROptimizerType::Pointer Roptimizer = ROptimizerType::New();
  RRegistrationType::Pointer Rregistration = RRegistrationType::New();

  Rregistration->SetMetric( Rmetric );
  Rregistration->SetOptimizer( Roptimizer );
  Rregistration->SetInterpolator( interpolator );

  Rregistration->SetTransform( Rtransform );
  Rregistration->SetFixedImage( fixedShrink->GetOutput() );
  Rregistration->SetMovingImage( movingShrink->GetOutput() );
  Rregistration->SetFixedImageRegion( fixedRegion );

  Rtransform->SetParameters( RfinalParameters );
  Rregistration->SetInitialTransformParameters( Rtransform->GetParameters() );

  OptimizerScalesType optimizerScales( Rtransform->GetNumberOfParameters() );

  const double translationScale = 1.0 / (2e4 *(it+1));
  const double scalingScale = 1.0 / 200;

  for(int i=0;i<9;i++){
    optimizerScales[i] = 1.0;
  }
  optimizerScales[0] = scalingScale;
  optimizerScales[4] = scalingScale;
  optimizerScales[8] = scalingScale;
  optimizerScales[9] = translationScale;
  optimizerScales[10] = translationScale;
  optimizerScales[11] = translationScale;
  Roptimizer->SetScales( optimizerScales );

  Roptimizer->SetMaximumStepLength( 0.5000  );
  Roptimizer->SetMinimumStepLength( 0.00001 );
  Roptimizer->SetNumberOfIterations( atoi(argv[5])*it );
  RObserver::Pointer Robserver = RObserver::New();
  Roptimizer->AddObserver( itk::IterationEvent(), Robserver );

  std::cout << std::endl << "Starting Rigid Registration" << std::endl;

  try {
    Rregistration->StartRegistration();}
  catch( itk::ExceptionObject & err ){
    std::cerr << "ExceptionObject caught during registration!" << std::endl;
    std::cerr << err << std::endl;
    return -1;}

  RfinalParameters = Rregistration->GetLastTransformParameters();

  const unsigned int numberOfIterations = Roptimizer->GetCurrentIteration();
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
  std::cout << "Center of rotation = " << std::endl << center << std::endl << std::endl;

	std::cout << "\n ----- " << it << "\\" << minlevel << "\n";
}


	//std::cout << "--------------- " << "TEST" << std::endl;
  RTransformType::Pointer finalTransform = RTransformType::New();
  finalTransform->SetCenter( Rtransform->GetCenter() );
  finalTransform->SetParameters( RfinalParameters );

	//std::cout << "--------------- " << "TEST 1" << std::endl;
  ResampleFilterType::Pointer Rresampler = ResampleFilterType::New();
  Rresampler->SetTransform( finalTransform );
  Rresampler->SetInterpolator( interpolator);
  Rresampler->SetInput( movingImageReader->GetOutput() );
  Rresampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  Rresampler->SetOutputOrigin(  fixedImage->GetOrigin() );
  Rresampler->SetOutputSpacing( fixedSpacing );
  Rresampler->SetDefaultPixelValue( 0 );
	//std::cout << "--------------- " << "TEST 1" << std::endl;

  try{
    Rresampler->Update();  }
  catch ( itk::ExceptionObject & err ){
    std::cerr << "An error occured while resampling the moving image after rigid registration" << std::endl;
    return EXIT_FAILURE;}


// Write the output image and transforms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//

		//std::cout << "--------------- " << "TEST 2" << std::endl;

#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
		itk::TransformFileWriterTemplate<float>::Pointer dWriter =
			itk::TransformFileWriterTemplate<float>::New();
		itk::TransformFileWriterTemplate<float>::Pointer iWriter =
			itk::TransformFileWriterTemplate<float>::New();
#else
		itk::TransformFileWriter::Pointer dWriter = itk::TransformFileWriter::New();
		itk::TransformFileWriter::Pointer iWriter = itk::TransformFileWriter::New();
#endif
		//std::cout << "--------------- " << "TEST 3" << std::endl;

  if(saveDTransform){  
	  // std::cout << "Exporting transformation parmeters to " << OutputTransformFileName << std::endl;
		//std::cout << finalTransform << std::endl;
		/*
	  RTransformType::MatrixType matrix = finalTransform->GetMatrix();
	  RTransformType::OffsetType offset = finalTransform->GetOffset();
      	  OutputTransform << "Translation" << std::endl << offset << std::endl << std::endl;
      	  OutputTransform << "Rotation" << std::endl << matrix;
		*/
		dWriter->SetInput(finalTransform);
		std::cout << "READY TO WRITE: " << std::endl;
		std::cout << argv[argcDTransform] << std::endl;
		dWriter->SetFileName(argv[argcDTransform]);
		try{
			dWriter->Update();
			std::cout << "Write direct transform as " << argv[argcDTransform] << std::endl;
		}catch( itk::ExceptionObject &err ){
			std::cerr << "Error while writing the direct transform as "
				<< argv[argcDTransform]
				<< std::endl;
		}
  }
	if( saveITransform ){
		std::cout << "READY TO WRITE: " << std::endl;
		std::cout << argv[argcITransform] << std::endl;
		iWriter->SetInput(finalTransform->GetInverseTransform());
		iWriter->SetFileName(argv[argcITransform]);
		try{
			iWriter->Update();
			std::cout << "Write inverse transform as " << argv[argcDTransform] << std::endl;
		}catch( itk::ExceptionObject &err ){
			std::cerr << "Error while writing the inverse transform as "
				<< argv[argcITransform]
				<< std::endl;
		}
	}

  IteratorType itr( fixedImage, fixedImage->GetLargestPossibleRegion() );
   itr.GoToBegin();
   while ( ! itr.IsAtEnd() ){
     const ImageType::IndexType index = itr.GetIndex();
     PixelType v = (Rresampler->GetOutput())->GetPixel(index);
     itr.Set( v );
     ++itr;
   }

  std::cerr << "Writing output image in "<< argv[3] << std::endl;    
    WriterType::Pointer outputWriter = WriterType::New();
    std::string outputname = argv[3];
    outputWriter->SetFileName( outputname.c_str() );
    outputWriter->SetInput( /*fixedImage*/ Rresampler->GetOutput() );
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


