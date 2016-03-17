/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  Perform a watershed segmentation
 *
 *        Version:  1.0
 *        Created:  09/02/15 10:25:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *


 * Watersheds segmentation
 * This technique is less sensitive
 * to user-defined thresholds than classic region-growing methods, and may be
 * better suited for fusing different types of features from different data
 * sets. The watersheds technique is also more flexible in that it does not 
 * produce a single image segmentation, but rather a hierarchy of 
 * segmentations from which a single region or set of regions can be 
 * extracted a-priori, using a threshold, or interactively, with the help of 
 * a graphical user interface
 *
 * The drawback of watershed segmentation is that it produces a region for
 * each local minimum in practice too many regions—and an over segmentation 
 * results. To alleviate this, we can establish a minimum watershed depth.
 *
 * Typically, the best results are obtained by preprocessing the original
 * image with an edge-preserving diffusion filter, such as one of the
 * anisotropic diffusion filters, or with the bilateral image filter.
 *
 * There are two parameters. Level controls watershed depth, and Threshold 
 * controls the lower thresholding of the input. Both parameters are set as
 * a percentage (0.0 - 1.0) of the maximum depth in the input image
 * =====================================================================================
 */


// my libraries
#include "../../dcmreader/dcmreader.h"
#include "../../dcmwriter/dcmwriter.h"
#include "../../console_tools/color.h"

#include "itkImage.h"
#include "itkWatershedImageFilter.h"
#include "itkCastImageFilter.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 2 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputImagePath"
			<<" [levels=0.1] [threshold=0.1] [PrincipalComponent=1]"
			<<std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argvInputImage = 1, argvOutputImage = 2,
							 argclevel = 3, argcThreshold = 4, argcPC=5;

	float level = 0.1, threshold = 0.1;
	bool usePC = true;
	if( argc > argclevel ){
		level = atof( argv[argclevel] );
	}
	if( argc > argcThreshold ){
		threshold = atof(argv[argcThreshold] );
	}
	if( argc > argcPC ){
		int tmp = atoi( argv[argcPC] );
		if( tmp == 1 )
			usePC = true;
		else
			usePC = false;
	}

	std::string msg;
	msg ="Watershed segmentation with parameters: level=";
	msg += SSTR(level);
	msg += " and threshold=" ;
	msg += SSTR(threshold);
	blueMessage( msg );

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
	typedef itk::WatershedImageFilter< ImageType >
																					WatershedFilterType;
	typedef itk::CastImageFilter< WatershedFilterType::OutputImageType,
																ImageType> CastImageFilterType;

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

	WatershedFilterType::Pointer watershed = WatershedFilterType::New();
	watershed->SetInput( image );
	watershed->SetLevel( level );
	watershed->SetThreshold( threshold );

	CastImageFilterType::Pointer cast = CastImageFilterType::New();
	cast->SetInput( watershed->GetOutput() );

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argvOutputImage] );
	writer.SetInput( cast->GetOutput() );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the volume as: "
			<< argv[argvOutputImage] << std::endl;
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}
	
	 msg = "Output volume wrote with ";
	msg += argv[argvOutputImage];
	finalMessage( msg );
	std::cout << std::endl;
}

