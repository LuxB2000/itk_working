/*
 * A complete pipeline that
 * 2- extract the brain
 * 3- measure the mean and std of the tissue along the slice
 */


// my libraries
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"

#include "itkImage.h"
#include "itkImageRegionIterator.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 3 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath"
			" PathToFileToWriteMeasurements"<<std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argcInputImage = 1,
							 argcPathToOutputFile = 2;

	// extract inputs
	std::string inputImagePath = argv[argcInputImage],
							outputFile = argv[argcPathToOutputFile];

	std::ofstream fileStream;
	fileStream.open(outputFile.c_str(), std::ofstream::app );
	if( !fileStream.is_open() ){
		std::cerr << "Impossible to write statistics in file: "
			<< outputFile << std::endl;
		exit(EXIT_FAILURE);
	}

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image, atlas;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;

	//
	// read the image
	//
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( inputImagePath.c_str() );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();
		image->DisconnectPipeline();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the input image with "
			<< inputImagePath << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}


	//
	// Measure the mean and std of the tissues
	//
	ImageType::SizeType sz = image->GetRequestedRegion().GetSize();
	// TODO: identify the tissues
	typedef std::vector<double> VectorType;
	VectorType mean = VectorType(sz[2],0), n = VectorType(sz[2],0), 
						 std = VectorType(sz[2],0);
	PixelType bkg = 100, pix;
	ImageType::IndexType ind;
	typedef itk::ImageRegionIteratorWithIndex< ImageType >
																													IteratorIndexType;
	IteratorIndexType it(image, image->GetRequestedRegion());

	// compute the mean
	for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
		// check if it's tissue
		pix = it.Get();
		if( pix > bkg ){
			ind = it.GetIndex();
			mean.at(ind[2]) += pix;
			n.at(ind[2]) += 1;
		}
	}

	for( int z=0; z<sz[2]; z++ ){
		mean.at(z) = mean.at(z) / n.at(z);
	}


	// compute the standard deviation
	for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
		// check if it's tissue
		pix = it.Get();
		if( pix > bkg ){
			ind = it.GetIndex();
			std.at(ind[2]) += pow( pix-mean.at(ind[2]), 2 );
		}
	}

	for( int z=0; z<sz[2]; z++ ){
		std.at(z) = sqrt( std.at(z) / n.at(z) );
	}

	//
	// Write results
	//
	fileStream << inputImagePath << "\n";
	for( int z=0; z<sz[2]; z++ ){
		fileStream << z << ", " << mean.at(z) <<", " << std.at(z) << "\n";
	}
	fileStream << "--------------------------------------------------------";
	fileStream << "--------------------------------------------------------\n";
}

