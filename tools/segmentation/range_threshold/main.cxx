/*
 * threshold an image, set OutsideValue to all values that are not in the 
 * range bound +/- (int) (range/2)
 *
 * 2015
 * Author: J Plumat
 *
 */

// my dcm reader
#include "./../dcmreader/dcmreader.h"
#include "./../dcmwriter/dcmwriter.h"

// general
#include "itkImage.h"
#include "itkThresholdImageFilter.h"

/*
 * MAIN
 */
int main( int argc, char *argv[] )
{
	std::cout<<" -- "<<argv[0]<<" -- "<<std::endl;
	bool zArtefact = false, yArtefact = false;
	// check parameters
  if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:" << argv[0]
			<< " InputImageFile"
			<< " outputImagefile bound range\n"
			<< "[outsidevalue, default=0]"
			<< std::endl;
    return EXIT_FAILURE;
    }

	int argcInputImageFile=1, argcOutputImageFile=2, 
			argcInputBound=3, argcInputRange=4, argcInputOutSideValue=5;

	// default value
	float outsideValue = 0.0;
	float bound = atof( argv[argcInputBound] );
	float range = atof( argv[argcInputRange] );

	if( argc == argcInputOutSideValue+1 ){
		outsideValue = atof( argv[argcInputOutSideValue] );
	}

	int lowerThresh = bound - range/2;
	int upperThresh = bound + range/2;
	std::cout << "Bound: " << bound << " range: " << range
		<< "lower thresh/upper thresh: " << lowerThresh << "/"
		<< upperThresh
		<< " outsideValue: " << outsideValue << std::endl;

	//
	// main type def
	//
  const unsigned int    Dimension = 3;
  typedef float					 OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension >  ImageType;
	typedef itk::ThresholdImageFilter <ImageType>
																			ThresholdImageFilterType;

	ImageType::Pointer image;


	//
	// Readers
	//
	typedef DcmReader<ImageType,Dimension> DCMReaderType;
	DCMReaderType reader =
													 DCMReaderType(argv[argcInputImageFile]);

	try{
		reader.Update();
		image = reader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr<<"Error while reading the input image with "
			<< argv[argcInputImageFile] <<std::endl;
		std::cerr<<err<<std::endl;
		exit(EXIT_FAILURE);
	}
	/*
	ImageType::Pointer output = ImageType::New();
	output->SetRegions( input->GetRequestedRegion() );
	output->SetSpacing( input->GetSpacing() );
	output->SetOrigin( input->GetOrigin() );
	output->Allocate();
	*/

	// thresholding
	ThresholdImageFilterType::Pointer thresholdFilter
																= ThresholdImageFilterType::New();
	thresholdFilter->SetInput( image );
	thresholdFilter->ThresholdOutside( lowerThresh, upperThresh );
	thresholdFilter->SetOutsideValue( outsideValue );


	//
	// Writer
	//
	typedef DcmWriter< ImageType, Dimension > MyDCMWriterType;
	MyDCMWriterType writer = MyDCMWriterType(  );
	writer.SetFileName( argv[argcOutputImageFile] );
	writer.SetInput( thresholdFilter->GetOutput() );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err ){
		std::cerr << "Error while writing the output image with "
			<< argv[argcOutputImageFile] << std::endl;
	}
	std::cout << "Output image wrote: " 
			<< argv[argcOutputImageFile] << std::endl;
	std::cout<<std::endl;

}//end main
