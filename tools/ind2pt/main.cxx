/*
 * change IMAGETYPE::POINTTYPE to IMAGETYPE::INDEXTYPE, given an image
 * input
 *	- image
 *	- pt1 pt2 pt3 TODO read values from a file, see IndexfileReader
 *	output is on console or TODO in file
 */


// my libraries
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"

#include "itkImage.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc != 5 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath ptx pty ptz"<<std::endl;
    return EXIT_FAILURE;
  }
	const unsigned int argvInputImage = 1, argcX = 2, argcY = 3, argcZ = 4;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	ImageType::IndexType ind;
	ImageType::PointType pt;

	// image types
	ImageType::Pointer image;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;

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

	ind[0] = atoi( argv[argcX] );
	ind[1] = atoi( argv[argcY] );
	ind[2] = atoi( argv[argcZ] );

	image->TransformIndexToPhysicalPoint( ind, pt );

std::cout << ind[0] << " " << ind[1] << " " << ind[2] << " --> " << pt[0] << " " << pt[1] << " " << pt[2] << std::endl;

	std::cout << std::endl;
}

