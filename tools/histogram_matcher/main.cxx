/*
 * transform image2 based on image1 histogram
 */


// my file readers
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"

// general
#include "itkImage.h"
#include "itkHistogramMatchingImageFilter.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 2 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImage1 InputImage2 outputImage2Modified"
			<<" [nbrOfPoints-in percent=20] [numberOfHistogramLevel=1024]"
			<<std::endl;
    return EXIT_FAILURE;
    }
	unsigned int argvInputImage1= 1 , argvInputImage2=2, argvOutputImage = 3,
							 argvNbrOfPoint=4, argvNbrOfHistLevel=5;
	unsigned int nbrOfHistLevel = 1024;
	float nbrOfPoints = 0.2;

	if( argc > argvNbrOfHistLevel ){
		nbrOfHistLevel = (unsigned int) atoi(argv[argvNbrOfHistLevel]);
	}
	if( argc > argvNbrOfPoint ){
		nbrOfPoints = atof( argv[argvNbrOfPoint] );
	}

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image1,image2;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;
	typedef itk::HistogramMatchingImageFilter<ImageType, ImageType
																	> MatchingFilterType;

	// read the image
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argvInputImage1] );
	try{
		dcmReader.Update();
		image1 = dcmReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the input image1 as: " <<
			argv[argvInputImage1]<< std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	MyDCMReaderType dcmReader2 = MyDCMReaderType();
	dcmReader2.SetFileName( argv[argvInputImage2] );
	try{
		dcmReader2.Update();
		image2 = dcmReader2.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the input image2 as: " <<
			argv[argvInputImage1]<< std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}


	ImageType::SizeType sz = image1->GetRequestedRegion().GetSize();
	double nbrOfPixels = sz[1]*sz[0]*sz[2];

	MatchingFilterType::Pointer matcher = MatchingFilterType::New();
	matcher->SetSourceImage( image2 );
	matcher->SetReferenceImage( image1 );
	matcher->SetNumberOfHistogramLevels( nbrOfHistLevel );
	matcher->SetNumberOfMatchPoints( (unsigned int) (nbrOfPoints*nbrOfPixels) );
	matcher->ThresholdAtMeanIntensityOn();

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argvOutputImage] );
	writer.SetInput( matcher->GetOutput() );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while the volume as: "
			<< argv[argvOutputImage] << std::endl;
	}
	
	std::cout << std::endl;
}

