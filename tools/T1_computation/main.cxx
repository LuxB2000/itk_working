/*
 * A simple example
 */


// my libraries
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"
#include <math.h>

#define PI 3.14159265

#include "itkImage.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 7 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " Image_FlipAngle1 Image_FlipeAngle2 OutputT1_Values";
		std::cerr << " flipAngle1 (in degree) flipAngle2 (in degree)";
		std::cerr << " TR (in ms)";
		std::cerr << std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argcImageFA1 = 1, argcImageFA2 = 2, argcOutputT1Values = 3,
							 argcFA1 = 4, argcFA2 = 5, argcTR = 6;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer imageFA1, imageFA2;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;

	// get inputs
	float a1 = atof(argv[argcFA1]), a2 = atof(argv[argcFA2]),
				TR = atof(argv[argcTR]);
	std::cout << "inputs : a1=" << a1 << " a2=" << a2 << " TR=" << TR << std::endl;

	// read the input images
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argcImageFA1] );
	try{
		dcmReader.Update();
		imageFA1 = dcmReader.GetOutput();// dcmReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}
	dcmReader.SetFileName( argv[argcImageFA2] );
	try{
		dcmReader.Update();
		imageFA2 = dcmReader.GetOutput();// dcmReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	// initiate the T1 volume
	ImageType::Pointer t1 = ImageType::New();
	t1->SetRegions( imageFA1->GetLargestPossibleRegion() );
	t1->SetSpacing( imageFA1->GetSpacing() );
	t1->SetOrigin( imageFA1->GetOrigin() );
	try{
		t1->Allocate();
	}catch( itk::ExceptionObject & err ){
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	// compute each T1 value
	ImageType::SizeType sz = imageFA1->GetLargestPossibleRegion().GetSize();
	typename ImageType::PointType pt;
	ImageType::PixelType  t1v = 0.0;
	double r =0.0,
				sa1 = sin(a1*PI/180), sa2 = sin(a2*PI/180),
				ca1 = cos(a1*PI/180), ca2 = cos(a2*PI/180);
	std::cout << "sa1: " << sa1 << " ca1: " << ca1 << " sa2: " << sa2 << " ca2: " << ca2 << std::endl;
	typename ImageType::IndexType ind, ind2, ind3;

	for( ind[2]=0; ind[2]<sz[2]; ind[2]++){
		for( ind[1]=0; ind[1]<sz[1]; ind[1]++){
			for( ind[0]=0; ind[0]<sz[0]; ind[0]++){

				imageFA1->TransformIndexToPhysicalPoint(ind,pt);
				imageFA2->TransformPhysicalPointToIndex(pt,ind2);
				t1->TransformPhysicalPointToIndex(pt,ind3);

				r = imageFA1->GetPixel(ind) / imageFA2->GetPixel(ind2);
				t1v = static_cast<ImageType::PixelType> 
					( TR / log((r*sa2*ca1-sa1*ca2)/(r*sa2-sa1)) );
					//( -TR / log((sa2-r*sa1)/(sa2*ca1-r*sa1*ca2)) );
				//std::cout << ind3 << " " << r << " " << t1v << std::endl;
				t1->SetPixel(ind3, t1v );
			}
		}
	}

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argcOutputT1Values] );
	writer.SetInput( t1 );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the volume as: "
			<< argv[argcOutputT1Values] << std::endl;
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}
	
	std::string msg = "Output volume wrote with ";
	msg += argv[argcOutputT1Values];
	finalMessage( msg );
	std::cout << std::endl;

}

