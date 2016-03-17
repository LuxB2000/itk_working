/*
 * A simple binary thresholding filter. Return a 0/1 binary map
 */


// my libraries
#include "../../dcmreader/dcmreader.h"
#include "../../dcmwriter/dcmwriter.h"
#include "../../console_tools/color.h"

#include "itkImage.h"
#include "itkBinaryThresholdImageFilter.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 4 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputImagePath"
			<<" upperThreshold [lowerThreshold, default=0]"
			<<std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argcInputImage = 1, argcOutputImage = 2, 
							 argcLowerThreshold = 4, argcUpperThreshold = 3;

	int upperThreshold = atoi(argv[argcUpperThreshold]),
			lowerThreshold = 0;

	if( argc > argcLowerThreshold ){
		lowerThreshold = atoi( argv[argcLowerThreshold] );
	}

	std::string msg = "Binary thresholding with lowerThreshold=" 
		+ SSTR(lowerThreshold)
		+ " upperThreshold=" + SSTR(upperThreshold);
	blueMessage(msg);

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;
	typedef itk::BinaryThresholdImageFilter<
				ImageType, ImageType > FilterType;

	// image types
	ImageType::Pointer image;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;

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

	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( image );
	filter->SetOutsideValue( 0 );
	filter->SetInsideValue( 1 );
	filter->SetLowerThreshold( lowerThreshold );
	filter->SetUpperThreshold( upperThreshold );

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argcOutputImage] );
	writer.SetInput( filter->GetOutput() );
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

