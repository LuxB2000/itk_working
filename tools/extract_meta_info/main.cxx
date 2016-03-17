/*
 * Extract the meta information relative to a TIFF file and write them
 * into a file.
 * Inputs:
 *	- TIff image path
 *	- output text file path
 *	- ID: an optional id that will be wrote on the beginning of the line,
 *			could be useful for exporting data to 'Excel'
 */


// my libraries
#include "../dcmreader/dcmreader.h"
#include "../console_tools/color.h"


int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 3 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputFilePath"
			<<" [ID< default='']"<<std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argvInputImage = 1, argcOutputFile = 2, argcID = 3;
	std::string inputFilePath = argv[argvInputImage],
							outputFilePath = argv[argcOutputFile],
							id = "",
							paramFilePath;

	if( argc > 3 ){
		id = argv[argcID];
	}

	std::ofstream fileStream;
	fileStream.open(outputFilePath.c_str(), std::ofstream::app );
	if( !fileStream.is_open() ){
		std::cerr << "Impossible to write statistics in file: "
			<< outputFilePath << std::endl;
		exit(EXIT_FAILURE);
	}

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;

	// read the image
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argvInputImage] );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();// dcmReader->GetOutput();
		paramFilePath = dcmReader.GetParamPath();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	if( paramFilePath.length() == 0 ){
		std::cerr << "We can't find the relative parameter file."
			<< " Have you set a DCM?" << std::endl;
		return( EXIT_FAILURE );
	}
	ReadTiffImageParameters parametersReader =
			ReadTiffImageParameters( paramFilePath );
			parametersReader.Read();
	TiffParameters *p = parametersReader.GetParamters();

	fileStream << id << ", " 
		<< p->time_complete.tm_year << "/" << p->time_complete.tm_mon << "/"
		<< p->time_complete.tm_mday 
		<< ", "
		<< p->time_complete.tm_hour << ":" << p->time_complete.tm_min << ":"
		<< p->time_complete.tm_sec
		<< ","
		<< p->tr
		<< ","
		<< p->te
		<< ","
		<< p->ti
		<< ","
		<< p->fliplist
		<< "\n";

	fileStream.close();

}

