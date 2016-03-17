/*
 * crop an image by giving inputs starting and ending slices in all 
 * directions
 * usage:
 * ./cropManually --startX firstSliceX --endX endSliceX 
 *                --startY firstSliceY --endY endSliceY
 *                --startZ firstSliceZ --endZ endSliceZ
 *                -i inputPath -o outputPath
 */


// my libraries
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"
#include "itkImage.h"
#include "itkRegionOfInterestImageFilter.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	
	const int defaultValue = 0;
	int startX = defaultValue, startY = defaultValue, startZ = defaultValue,
			endX=defaultValue, endY=defaultValue, endZ=defaultValue,
			c = 0;
	std::string msg="",defaultStr="";
	std::string inputPath=defaultStr, outputPath=defaultStr;

	std::cout << argc << std::endl;


	msg="inputs: \n";
	
	// get inputs
	static struct option long_options[] =
			{
				// x options
				{"startX", required_argument,  0, 'a'},
				{"endX",   required_argument,  0, 'b'},
				// y options
				{"startY", required_argument,  0, 'c'},
				{"endY",   required_argument,  0, 'd'},
				// z options
				{"startZ", required_argument,  0, 'e'},
				{"endZ",   required_argument,  0, 'f'},
				// input and output volumes
				{"input",   required_argument, 0, 'i'},
				{"output",  required_argument, 0, 'o'},
				{0,                 0        , 0,  0 }
			};
	/* getopt_long stores the option index here. */
	int option_index = 0;

	while(1){
		c = getopt_long (argc, argv, "a:b:c:d:e:f:i:o:",
										 long_options, &option_index);

		// no more options, break the loop
		if( c== -1 ) break;

		switch(c){
			case 'a' :
				startX = atoi(optarg);
				msg += "\t-startX " + SSTR(startX) + "\n";
				break;
			case 'b' :
				endX = atoi(optarg);
				msg += "\t-endX " + SSTR(endX) + "\n";
				break;
			case 'c' :
				startY = atoi(optarg);
				msg += "\t-startY " + SSTR(startY) + "\n";
				break;
			case 'd' :
				endY = atoi(optarg);
				msg += "\t-endY " + SSTR(endY) + "\n";
				break;
			case 'e' :
				startZ = atoi(optarg);
				msg += "\t-startZ " + SSTR(startZ) + "\n";
				break;
			case 'f' :
				endZ = atoi(optarg);
				msg += "\t-endZ " + SSTR(endZ) + "\n";
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

	// print inputs message
	blueMessage(msg);

	// Check input parameters
	if( strcmp(inputPath.c_str(),defaultStr.c_str())==0 ||
			strcmp(outputPath.c_str(),defaultStr.c_str())==0 )
	{
		error("Missing Parameters");
		std::cerr << "Usage: -i InputPath -o OutputPath";
		std::cerr << "--startX x1 --endX x2 --startY y1 --endY y2";
		std::cerr << "--startZ z1 --endZ z2";
		std::cerr << std::endl;
		exit( EXIT_FAILURE);
	}

	if( startX<0 || startY<0 || startZ<0 || endX<0 || endY<0 || endZ<0 ){
		error("Slice number can't be negative");
		exit(EXIT_FAILURE);
	}


	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;
	typedef ImageType::SizeType SizeType;

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

	SizeType inputSize = image->GetLargestPossibleRegion().GetSize();

	// check if endSlice values are correct
	if(endX>=inputSize[0] || endX==defaultValue){
		endX = inputSize[0];
	}
	if(endY>=inputSize[1] || endY==defaultValue){
		endY = inputSize[1];
	}
	if(endZ>=inputSize[2] || endZ==defaultValue){
		endZ = inputSize[2];
	}

	if( endX<=startX || endY<=startY || endZ<=startZ ){
		error("Start slice can not be higher than end slices.");
		exit(EXIT_FAILURE);
	}

	typedef itk::RegionOfInterestImageFilter< ImageType,
			ImageType > ExtractorFilter;
	ExtractorFilter::Pointer extractor = ExtractorFilter::New();
	SizeType outputSize;
	outputSize[0] = endX-startX;
	outputSize[1] = endY-startY;
	outputSize[2] = endZ-startZ;
	ImageType::IndexType start;
	start[0] = startX;
	start[1] = startY;
	start[2] = startZ;

	std::cout << inputSize << std::endl;
	msg ="output size: [" + SSTR(outputSize[0]) + ", " + SSTR(outputSize[1])
		+ "," + SSTR(outputSize[2]) + "] start: [" + 
		SSTR(start[0]) + "," + SSTR(start[1]) + "," + SSTR(start[2]) + "]" ;
	blueMessage(msg);

	ImageType::RegionType outputRegion;
	outputRegion.SetSize(outputSize);
	outputRegion.SetIndex(start);

	extractor->SetInput( image );
	extractor->SetRegionOfInterest( outputRegion );

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( outputPath.c_str() );
	writer.SetInput( extractor->GetOutput() );
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

