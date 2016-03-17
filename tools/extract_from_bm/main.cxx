/*
 * Extract part of image based on a Binary Mask (BM).
 * The BM may contains different values
 * The program can extract positions specified by BM or extract positions 
 * not specify by BM.
 *
 * example:
 * ./extractBM volumeToExtract.mha BMSpecifyPositions.mha volume_extracted.mha
 * ./extractBM volumeToExtract.mha BMSpecifyPositions.mha volume_extracted.mha 1 0 294
 */


// my file readers
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
// displaying results
#include "../console_tools/color.h"

// general
#include "itkImage.h"
#include "itkImageRegionIterator.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
	unsigned int argvInputImage = 1, argvInputBM=2, argvOutputImage = 3,
							 argvKeepPos=4, argvBMPixVal=6, argvOutPutPix=5;

  if( argc < (argvOutputImage+1) )
    {
			error( "Missing Parameters ");
			std::cerr << "Usage:" << argv[0];
			std::cerr << " InputImagePath binaryMaskPath outputImagePath\n"
			<< "[Keep position specified by BM 1/0, default=1]\n"
			<< "[OutputPixel, default=0]\n"
			<< "[BinaryMaskPixel, default only one pixel is present]\n"
			<<std::endl;
    return EXIT_FAILURE;
    }

	//
	// Check inputs
	//
	bool extractPosBM = true, oneLabelInBM=true;
	float defautlOutputPix = 0.0, labelInBM=0.0;
	if( argc >= (argvBMPixVal+1) ){
		labelInBM = atof( argv[argvBMPixVal] );
		oneLabelInBM = false;
	}
	if( argc >= (argvOutPutPix + 1) ){
		defautlOutputPix = atof( argv[argvOutPutPix] );
	}
	if( argc >= (argvKeepPos+1) ){
		int b = (atoi(argv[argvKeepPos]));
		if( b!=0 && b!=1){
			std::cerr << "ERROR argv["<<argvKeepPos<<"] should be 0 or 1" <<
				std::endl;
			return( EXIT_FAILURE );
		}
		extractPosBM = (b==1) ? true : false;
	}

	std::cout << "Inputs: keep pos specified by BM=" << 
		extractPosBM ;
	std::cout << "- default output pixel=" << defautlOutputPix 
		<< " - chosen label=";
	if( oneLabelInBM )
		std::cout << "only one label should be present ";
	else
		std::cout << labelInBM << " ";
	std::cout << std::endl;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image, mask;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;
	typedef itk::ImageRegionIterator<ImageType>			 IteratorType;


	// loops variables
	PixelType p;
	unsigned int i=0, s=0;

	// read the image
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argvInputImage] );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the input volume as: "
			<< argv[argvInputImage] << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	// read the binary mask
	MyDCMReaderType maskReader = MyDCMReaderType();
	maskReader.SetFileName( argv[argvInputBM] );
	try{
		maskReader.Update();
		mask = maskReader.GetOutput();
	}catch( itk::ExceptionObject &err ){
		std::cerr << "Error while reading the input Binary Mask as: "
			<< argv[argvInputBM] << std::endl;
	}


	// parse the mask to find all pixel values
	IteratorType itM( mask, mask->GetRequestedRegion() );
	std::vector< PixelType > bmPixel = std::vector<PixelType>();
	for( itM.GoToBegin(); !itM.IsAtEnd(); ++itM ){
		p = itM.Get();
		if( p != 0 ){
			// check if the value is not already present in the vector
			if(std::find(bmPixel.begin(), bmPixel.end(), p) == bmPixel.end()) {
				bmPixel.push_back( p );
			}
		}
	}
	s = bmPixel.size();

	// check if BM is not empty
	if( s==0 ){
		std::cerr << "ERROR the Binary mask is empty !" << std::endl;
		return( EXIT_FAILURE );
	}

	if( oneLabelInBM && s==1 ){
		labelInBM = (PixelType) bmPixel.at(0);
	}

	if( oneLabelInBM && s>1 ){
		// we expect one label but more are found
		std::cerr << "ERROR: no label is specified in BM and we found more: ";
		for( i=0; i<bmPixel.size(); i++ ){
			std::cerr << bmPixel.at(i);
			if( i!=(s-1) ) std::cerr << ", ";
		}
		std::cerr << std::endl;
		std::cerr << "I don't know what to do. Please specify which label must"
			<< " be chosen."<< std::endl;
	}

	// check if the input label is present in BM
	if( !oneLabelInBM )
		if(std::find(bmPixel.begin(), bmPixel.end(), labelInBM)
				== bmPixel.end()) {
			std::cerr << "ERROR the specified label (" << labelInBM <<
				") is not in present in the BM.\n"
				<<"We found: ";
			for( i=0; i<bmPixel.size(); i++ ){
				std::cerr << bmPixel.at(i);
				if( i!=(s-1) ) std::cerr << ", ";
			}
			std::cerr << std::endl;
			return( EXIT_FAILURE );
		}


	//
	// Extraction
	//
	IteratorType itV( image, image->GetRequestedRegion() );
	p = (PixelType) defautlOutputPix;
	for( itM.GoToBegin(), itV.GoToBegin(); !itM.IsAtEnd(), !itV.IsAtEnd();
			++itV, ++itM ){
		if( itM.Get() == labelInBM ){
			if( !extractPosBM ){
				itV.Set( p );
			}
		}else{
			if( extractPosBM ){
				itV.Set( p );
			}
		}
	}

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argvOutputImage] );
	writer.SetInput( image );
	try{
		writer.Update();
		finalMessage( "Output Image wrote as: " +
				std::string(argv[argvOutputImage]) );
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the volume as: "
			<< argv[argvOutputImage] << std::endl;
		std::cerr << err << std::endl;
	}
	
	std::cout << std::endl;
}

