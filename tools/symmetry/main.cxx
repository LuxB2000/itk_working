/*
 * A reflexive symmetry applied on a certain axis. The axis is place on the
 * center of the image.
 */


// my file readers
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"

// general
#include "itkImage.h"
#include <itkImageConstIteratorWithIndex.h>
#include "itkImageRegionIteratorWithIndex.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 2 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputImagePath [axis]"<<std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argvInputImage = 1, argvOutputImage = 2, argcAxis = 3;
	int axis = 1; // 0: X axis, 1: Y axis, 2: Z axis
	if( argc >= argcAxis ){
		axis = atoi( argv[argcAxis] );
		if( axis!= 0 && axis!=1 && axis != 2){
			error("Wrong axis");
			std::cerr << "Axis should be 0 (X axis), 1 (Y axis) or 2 (Z axis)" 
				<< " input=" << axis << std::endl;
			exit( EXIT_FAILURE );
		}
	}
	std::cout << "Input Axis: " << axis << std::endl;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;
	typedef itk::ImageRegionIteratorWithIndex< ImageType > Iterator;

	// image types
	ImageType::Pointer image;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;

	// read the image
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argvInputImage] );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();// dcmReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	ImageType::Pointer output = ImageType::New();
	output->SetRegions( image->GetRequestedRegion() );
	output->SetOrigin( image->GetOrigin() );
	output->SetSpacing( image->GetSpacing() );
	try{
		output->Allocate();
	}catch( itk::ExceptionObject &err ){
		error("Error while allocating the output image.");
		std::cerr << err << std::endl;
	}

	Iterator outputIt(output, output->GetRequestedRegion() );
	ImageType::IndexType requestedIndex =
												image->GetRequestedRegion().GetIndex();
	ImageType::SizeType requestedSize =
												image->GetRequestedRegion().GetSize();
	for ( outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++outputIt)
	{
		ImageType::IndexType idx = outputIt.GetIndex();
		idx[axis] = requestedIndex[axis] + requestedSize[axis] - 1 - idx[axis];
		outputIt.Set( image->GetPixel(idx) );
	}



	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argvOutputImage] );
	writer.SetInput( output );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::string msg( "Error while writing the volume as: ");
			msg += argv[argvOutputImage];
		error(msg);
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}
	
	std::string msg = "Output volume wrote with ";
	msg += argv[argvOutputImage];
	finalMessage( msg );
	std::cout << std::endl;
}

