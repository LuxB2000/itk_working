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
  if( argc < 4 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputImagePath";
		std::cerr << " ReferenceImagePath";
    return EXIT_FAILURE;
  }
	unsigned int argcInputImage = 1, argcOutputImage = 2,
							 argcRefImagePath = 3;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image, imageRef;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;

	// read the images
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argcInputImage] );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	dcmReader.SetFileName( argv[argcRefImagePath] );
	try{
		dcmReader.Update();
		imageRef = dcmReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	ImageType::SpacingType refSpacing   = imageRef->GetSpacing();
	ImageType::SpacingType inputSpacing = image->GetSpacing();

	ImageType::SizeType outputSize, refSize  = imageRef->GetLargestPossibleRegion().GetSize();
	ImageType::SizeType inputSize= image->GetLargestPossibleRegion().GetSize();

	std::string msg = "";
	msg += "\ninput size " + SSTR(inputSize[0]) + " "
		+ SSTR(inputSize[1]) + " " + SSTR(inputSize[2]);
	msg = "\ninput spacing " + SSTR(inputSpacing[0]) + " "
		+ SSTR(inputSpacing[1]) + " " + SSTR(inputSpacing[2]);
	msg += "\nref spacing " + SSTR(refSpacing[0]) + " "
		 + SSTR(refSpacing[1]) + " " + SSTR(refSpacing[2]);

	// outputSpacing[i] = inputSpacing[i] * (static_cast<double>(inputSize[i]) / static_cast<double>(outputSize[i]));
	for(unsigned int i=0; i<3; i++)
		outputSize[i] = static_cast<int>( (inputSpacing[i]/refSpacing[i]) * inputSize[i] );

	msg += "\noutput size " + SSTR(outputSize[0]) + " " 
	 + SSTR(outputSize[1]) + " "  + SSTR(outputSize[2]);

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
  resample->SetOutputSpacing(refSpacing);
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

