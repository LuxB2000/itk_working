/*
 * Read a SHORT meta image and write it as Ushort
 */


#include "itkImageFileReader.h"
// my file readers
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"

// general
#include "itkImage.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 2 )
    {
    error( "Missing Parameters" );
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputImagePath"<<std::endl;
    return EXIT_FAILURE;
    }
	unsigned int argvInputImage = 1, argvOutputImage = 2;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
	typedef short		InputPixelType;
  typedef  float          PixelType;
	typedef itk::Image<InputPixelType,Dimension> InputImageType;

	// image types
	InputImageType::Pointer image;
	typedef DcmWriter<InputImageType, Dimension>    MyDCMWriterType;

	// read the image
	typedef itk::ImageFileReader<InputImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( argv[argvInputImage] );
	try{
		reader->Update();
		image = reader->GetOutput();// dcmReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argvOutputImage] );
	writer.SetInput( image );
	try{
		writer.Update();
		finalMessage( "Output saved as: " + std::string(argv[argvOutputImage]));
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while the volume as: "
			<< argv[argvOutputImage] << std::endl;
	}
	
	std::cout << std::endl;
}

