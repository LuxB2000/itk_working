/*
 * A simple example that uses getopt
 *
 * from https://www.gnu.org/software/libc/manual/html_node/Using-Getopt.html#Using-Getopt
 * int getopt (int argc, char *const *argv, const char *options)
 * The options argument is a string that specifies the option characters that are valid for this program. An option character in this string can be followed by a colon (‘:’) to indicate that it takes a required argument. If an option character is followed by two colons (‘::’), its argument is optional; this is a GNU extension.
 *
 * Example usage
 * ./readAndWriteGetopt  -i /some/path -o /other/path -a
 * ./readAndWriteGetopt  -i /some/path -o /other/path --someFlag
 */


// my libraries
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"
#include "itkImage.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	
	int aFlag = 0, c = 0;
	std::string msg="",defaultStr="";
	std::string inputPath=defaultStr, outputPath=defaultStr;

	msg="inputs: \n";
	
	// get inputs
	static struct option long_options[] =
	{
		/* These options don’t set a flag.
			 We distinguish them by their indices. */
		{"someFlag",     no_argument,  0, 'a'},
		// input and output volumes
		{"input",  required_argument,  0, 'i'},
		{"output",  required_argument, 0, 'o'},
		{0,                 0        , 0,  0 }
	};
	/* getopt_long stores the option index here. */
	int option_index = 0;

	while(1){
		c = getopt_long (argc, argv, "i:o:a",
										 long_options, &option_index);

		// no more options, break the loop
		if( c== -1 ) break;

		switch(c){
			case 'a' :
				aFlag = 1;
				msg += "\t-a " + SSTR(aFlag) + "\n";
				break;
			// input volume
			case 'i' :
				inputPath = std::string(optarg);
				msg += "\t-i " + inputPath + "\n";
				break;
			// output volume
			case 'o' :
				outputPath = std::string(optarg);
				msg += "\t-o " + outputPath + "\n";
				break;
			default :
				abort();
		}
	}

	if( strcmp(inputPath.c_str(),defaultStr.c_str())==0 ||
			strcmp(outputPath.c_str(),defaultStr.c_str())==0 )
	{
		error("Missing Parameters");
		std::cerr << "Usage: -i InputPath -o OutputPath";
		std::cerr << std::endl;
		exit( EXIT_FAILURE);
	}


	blueMessage(msg);

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

	// read the image
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( inputPath.c_str() );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();// dcmReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( outputPath.c_str() );
	writer.SetInput( image );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the volume as: "
			<< outputPath << std::endl;
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}
	
	msg = "Output volume wrote with ";
	msg += outputPath;
	finalMessage( msg );
	std::cout << std::endl;
}

