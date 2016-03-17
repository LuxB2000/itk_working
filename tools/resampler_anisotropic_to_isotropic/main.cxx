/*
 * A simple example
 */


// my libraries
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"

#include "itkImage.h"
#include "itkResampleImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkIntensityWindowingImageFilter.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 3 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputAnisotropicImagePath outputIsotropicImagePath";
		std::cerr << std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argcInputImage = 1, argcOutputImage = 2;

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
	typedef itk::RecursiveGaussianImageFilter<
										ImageType, ImageType > GaussianFilterType;
	GaussianFilterType::Pointer smootherX = GaussianFilterType::New();
	GaussianFilterType::Pointer smootherY = GaussianFilterType::New();
	typedef itk::ResampleImageFilter< ImageType, 
														ImageType > ResampleFilterType;
	ResampleFilterType::Pointer resamplerFilter = ResampleFilterType::New();
	typedef itk::IdentityTransform < double, Dimension > TransformType;
	TransformType::Pointer transform = TransformType::New();
	typedef itk::LinearInterpolateImageFunction< ImageType, double >
															InterpolateType;
	InterpolateType::Pointer interpolator = InterpolateType::New();

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

	// get the input (anisotropic) spacing
	const ImageType::SpacingType& inputSpacing = image->GetSpacing();
	const double isoSpacing = std::sqrt( inputSpacing[2] * inputSpacing[0] );
	ImageType::SizeType outputSize, inputSize = image->GetLargestPossibleRegion().GetSize();
	typedef ImageType::SizeType::SizeValueType SizeValueType;
	const double dx = inputSize[0] * inputSpacing[0] / isoSpacing;
	const double dy = inputSize[1] * inputSpacing[1] / isoSpacing;
	const double dz = ( inputSize[2] - 1 ) * inputSpacing[2] / isoSpacing;

	ImageType::SpacingType outSpacing;
	outSpacing[0] = isoSpacing;
	outSpacing[1] = isoSpacing;
	outSpacing[2] = isoSpacing;

	outputSize[0] = static_cast<SizeValueType>(dx);
	outputSize[1] = static_cast<SizeValueType>(dy);
	outputSize[2] = static_cast<SizeValueType>(dz);


	// Set smoother
	smootherX->SetSigma( isoSpacing );
	smootherX->SetDirection( 0 );
	smootherX->SetInput( image );

	smootherY->SetSigma( isoSpacing );
	smootherY->SetDirection( 1 );
	smootherY->SetInput( smootherX->GetOutput() );

	// resampling
	transform->SetIdentity();
	resamplerFilter->SetInput( image );
	resamplerFilter->SetTransform( transform );
	resamplerFilter->SetInterpolator( interpolator );
	resamplerFilter->SetOutputSpacing( outSpacing );
	resamplerFilter->SetOutputOrigin( image->GetOrigin() );
	resamplerFilter->SetOutputDirection( image->GetDirection() );
	resamplerFilter->SetSize( outputSize );
	resamplerFilter->Update();

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argcOutputImage] );
	writer.SetInput( resamplerFilter->GetOutput() );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the volume as: "
			<< argv[argcOutputImage] << std::endl;
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}
	
	std::string msg = "Output volume wrote with ";
	msg += argv[argcOutputImage];
	finalMessage( msg );
	std::cout << std::endl;
}

