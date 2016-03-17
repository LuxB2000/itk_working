/*
 * A simple example
 */


// my libraries
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"

#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkImage.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 6 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputImagePath";
		std::cerr << " outputXSize outputYSize outputZSize";
    return EXIT_FAILURE;
  }
	unsigned int argcInputImage = 1, argcOutputImage = 2,
							 argcXSize = 3, argcYSize = 4, argcZSize = 5;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;

	// read user inputs
	ImageType::SizeType outputSize;
	outputSize.Fill(0);
	outputSize[0] = atoi( argv[argcXSize] );
	outputSize[1] = atoi( argv[argcYSize] );
	outputSize[2] = atoi( argv[argcZSize] );

	// read the image
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argcInputImage] );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();// dcmReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	ImageType::SpacingType outputSpacing,inputSpacing = image->GetSpacing();
	ImageType::SizeType inputSize = image->GetLargestPossibleRegion().GetSize();
	std::cout << "Test : " << image->GetLargestPossibleRegion().GetSize() << std::endl;
	std::string msg = "";
	msg += "\ninput size " + SSTR(inputSize[0]) + " "
		+ SSTR(inputSize[1]) + " " + SSTR(inputSize[2]);
	msg = "\ninput spacing " + SSTR(inputSpacing[0]) + " "
		+ SSTR(inputSpacing[1]) + " " + SSTR(inputSpacing[2]);
	msg += "\noutput size " + SSTR(outputSize[0]) + " "
		 + SSTR(outputSize[1]) + " " + SSTR(outputSize[2]);

	for(unsigned int i=0; i<3; i++)
	outputSpacing[i] = inputSpacing[i] * (static_cast<double>(inputSize[i]) / static_cast<double>(outputSize[i]));

	msg += "\noutput spacing " + SSTR(outputSpacing[0]) + " " 
	 + SSTR(outputSpacing[1]) + " "  + SSTR(outputSpacing[2]);

	blueMessage(msg);

	typedef itk::NearestNeighborInterpolateImageFunction<
                       ImageType, double >  InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	// resampling
	typedef itk::IdentityTransform<double, 3> TransformType;
  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleImageFilterType;
	ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();
  resample->SetInput(image);
	resample->SetInterpolator( interpolator );
  resample->SetSize(outputSize);
  resample->SetOutputSpacing(outputSpacing);
	resample->SetOutputOrigin( image->GetOrigin() );
	resample->SetOutputDirection( image->GetDirection() );
  resample->SetTransform(TransformType::New());
  resample->UpdateLargestPossibleRegion();

	/*
	ImageType::IndexType ind; ind[0] = 10; ind[1] = 10; ind[2] = 4;
	std::cout << "Test 2.1: " << image->GetPixel(ind) << std::endl;
	std::cout << "Test 2 : " << resample->GetOutput()->GetPixel(ind) << std::endl;
	*/
	
	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argcOutputImage] );
	writer.SetInput( resample->GetOutput() );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the volume as: "
			<< argv[argcOutputImage] << std::endl;
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}
	
	msg = "Output volume wrote with ";
	msg += argv[argcOutputImage];
	finalMessage( msg );
	std::cout << std::endl;
}

