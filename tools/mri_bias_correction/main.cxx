/*
 * Correct the MRI bias
 * 
 * INPUTS
 * ------
 * * Input image path
 * * Output image path
 * OUTPUTS
 * -------
 * * If no error crop version is written on HD
 *
 */


// my file readers
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"

// general
#include "itkImage.h"

#include "itkGaussianMixtureModelComponent.h"
#include "itkExpectationMaximizationMixtureModelEstimator.h"
#include "itkListSample.h"

// MRI bias correction
#include <itkMRIBiasFieldCorrectionFilter.h>

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 2 )
    {
    error( "Missing Parameters " );
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputCorrImagePath "<<std::endl;
    return EXIT_FAILURE;
    }
	unsigned int argvInputImage = 1, argvOutputImage = 2;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	std::string msg;

	// image types
	typedef itk::Image< PixelType, Dimension > ImageType;
	typedef itk::Image< bool, Dimension >      MaskType;
	ImageType::Pointer image, corr_image;
	ImageType::Pointer cropImage;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension >   MyDCMWriterType;
	typedef	itk::MRIBiasFieldCorrectionFilter< ImageType, ImageType, ImageType>
																						MRIBiasFieldCorrectionFilterType;

	//
	// read the mask image with my reader
	//
	try{
		MyDCMReaderType reader = MyDCMReaderType( argv[argvInputImage] );
		reader.Update();
		image = reader.GetOutput();
		image->DisconnectPipeline();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading input volume as " <<
			argv[argvInputImage]  << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	//
	// MRI Bias correction
	//
	typedef itk::Vector< double, 1 > MeasurementVectorType;
	typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;

	// prepare the mask on the whole image
	ImageType::RegionType region = image->GetLargestPossibleRegion();
	ImageType::Pointer mask = ImageType::New();
	mask->SetRegions(region);
	mask->Allocate();
	
	itk::ImageRegionIterator<ImageType> imageIterator(mask,region);
	while(!imageIterator.IsAtEnd())
  {
		imageIterator.Set(1);
		++imageIterator;
	}

	// statistics for the two classes: foreground and background
	// bkg is assume to have values below 1000, while foreground is higher

	unsigned int numberOfClasses = 2;
	typedef itk::Array< double > ParametersType;
	ParametersType param(numberOfClasses),	mean(numberOfClasses),
																					std(numberOfClasses);
	std::vector< ParametersType > initialParameters( numberOfClasses );
	typedef itk::Statistics::ExpectationMaximizationMixtureModelEstimator<
		SampleType > EstimatorType;
	typedef itk::Statistics::GaussianMixtureModelComponent< SampleType >
																											ComponentType;

	SampleType::Pointer sample = SampleType::New();
	sample->SetMeasurementVectorSize( 1 ); // length of measurement vectors
																				//in the sample, here pixel intensity
	itk::ImageRegionIterator<ImageType> imageIterator2(image,region);
	while(!imageIterator2.IsAtEnd())
  {
		sample->PushBack(imageIterator2.Get());
		++imageIterator2;
	}
	param[0] = 600; // initial bkg values
	param[1] = 800;
	initialParameters[0] = param;
	param[0] = 5000; // initial foreground values
	param[1] = 800;
	initialParameters[1] = param;
	std::vector< ComponentType::Pointer > components;
	for ( unsigned int i = 0; i < numberOfClasses; i++ )
	{
		components.push_back( ComponentType::New() );
		(components[i])->SetSample( sample );
		(components[i])->SetParameters( initialParameters[i] );
	}


	EstimatorType::Pointer estimator = EstimatorType::New();
	estimator->SetSample( sample );
	estimator->SetMaximumIteration( 200 );

	itk::Array< double > initialProportions(numberOfClasses);
	initialProportions[0] = 0.5;
	initialProportions[1] = 0.5;

	estimator->SetInitialProportions( initialProportions );
	for ( unsigned int i = 0; i < numberOfClasses; i++)
	{
		estimator->AddComponent( (ComponentType::Superclass*)
		(components[i]).GetPointer() );
	}
	try{
		estimator->Update();
	}catch( itk::ExceptionObject &err ){
		std::cout << "Error while estimating the mean and std of classes." 
			<< std::endl;
		std::cout << err << std::endl;
		return (EXIT_FAILURE);
	}

	// print EM estimations
	msg = "";
	for ( unsigned int i = 0; i < numberOfClasses; i++ )
	{
	msg+= std::string("Cluster[") + SSTR( i ) + std::string( "]\n");
	msg+= std::string("Parameters:")+ 
				SSTR((components[i])->GetFullParameters()) + std::string("\n");
	msg+= std::string("Proportion: ") + 
				SSTR(estimator->GetProportions()[i]) + std::string("\n\n");
	}
	blueMessage( msg );

	// extract information for corrective filter
	mean[0] = components[0]->GetFullParameters()[0];
	std[0] = components[0]->GetFullParameters()[1];
	mean[1] = components[1]->GetFullParameters()[0];
	std[1] = components[1]->GetFullParameters()[1];

	// Correction filter
	MRIBiasFieldCorrectionFilterType::Pointer filter =
																	MRIBiasFieldCorrectionFilterType::New();
	filter->SetUsingSlabIdentification( false ); // default: false
	filter->SetUsingInterSliceIntensityCorrection( true); // default: true
	filter->SetUsingBiasFieldCorrection( true ); // default: true
	filter->SetNumberOfLevels( 2 ); // default: 2

	filter->SetInput( image );
	filter->SetInputMask( mask );
	filter->SetTissueClassStatistics( mean, std );
	
	/*
	typedef itk::CommandIterationUpdate<MRIBiasFieldCorrectionFilterType>
		CommandType;
  CommandType::Pointer observer = CommandType::New();
  filter->AddObserver( itk::IterationEvent(), observer );
	*/
	filter->Update();

	//
	// write corrected volume
	//
	try{
		MyDCMWriterType writer = MyDCMWriterType( argv[argvOutputImage] );
		writer.SetInput( filter->GetOutput() );
		//writer.SpecifyOutputType( 2 ); // write as dicom stack
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing output volume as "
			<< argv[argvOutputImage] << std::endl;
		std::cerr << err << std::endl;
	}
	msg = "Corrected bias image wrote with success to " + std::string( argv[argvOutputImage]);
	finalMessage( msg );
	std::cout << std::endl;

}

