/*
 * Apply a gaussian smoothing to an image
 *
 * 2015
 * Author: J Plumat
 *
 */

// my dcm reader
#include "./../dcmreader/dcmreader.h"
#include "./../dcmwriter/dcmwriter.h"

// general
#include "itkImage.h"
#include "itkGaussianOperator.h"

// iterators
#include "itkImageRegionIterator.h"
//#include "itkImageRegionConstIterator.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkRescaleIntensityImageFilter.h"


/*
 * MAIN
 */
int main( int argc, char *argv[] )
{
	std::cout<<" -- "<<argv[0]<<" -- "<<std::endl;
	bool zArtefact = false, yArtefact = false;
	// check parameters
  if( argc < 2 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:" << argv[0]
			<< " InputImageFile"
			<< " outputImagefile"
			<<" [standardDeviation]"
			<< std::endl;
    return EXIT_FAILURE;
    }

	int argcInputImageFile=1, argcOutputImageFile=2, 
			argcInputStandardDeviation=3;

	// default value
	float stdDev = 2.5;

	if( argc == argcInputStandardDeviation+1 ){
		stdDev = atof( argv[argcInputStandardDeviation] );
	}

	std::cout << "Standard deviation: " << stdDev << std::endl;

	//
	// main type def
	//
  const unsigned int    Dimension = 3;
  typedef float					 OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension >  ImageType;
	typedef itk::ImageRegionIterator< ImageType > IteratorType;
	typedef itk::GaussianOperator< OutputPixelType, Dimension >
																										GaussianOperatorType;
	 typedef itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< 
																						ImageType > FaceCalculatorType;
	itk::NeighborhoodInnerProduct<ImageType> innerProduct;

	FaceCalculatorType faceCalculator;
	FaceCalculatorType::FaceListType faceList;
	FaceCalculatorType::FaceListType::iterator fit;

	typedef itk::ConstNeighborhoodIterator< ImageType > 
																	NeighborhoodIteratorType;
	typedef itk::ImageRegionIterator< ImageType> IteratorType;
	IteratorType out;
	NeighborhoodIteratorType it;


	ImageType::Pointer input;


	//
	// Readers
	//
	typedef DcmReader<ImageType,Dimension> DCMReaderType;
	DCMReaderType reader =
													 DCMReaderType(argv[argcInputImageFile]);

	try{
		reader.Update();
		input = reader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr<<"Error while reading the input image with "
			<< argv[argcInputImageFile] <<std::endl;
		std::cerr<<err<<std::endl;
		exit(EXIT_FAILURE);
	}
	ImageType::Pointer output = ImageType::New();
	output->SetRegions( input->GetRequestedRegion() );
	output->SetSpacing( input->GetSpacing() );
	output->SetOrigin( input->GetOrigin() );
	output->Allocate();


	//
	// smoothing
	//
	GaussianOperatorType gaussianOperator;
	gaussianOperator.SetVariance( stdDev * stdDev );
	for( unsigned int d=0; d<ImageType::ImageDimension; d++ ){
		std::cout << "smooth dimension: " << d << std::endl;
		gaussianOperator.SetDirection(d);
		gaussianOperator.CreateDirectional();
		 faceList = faceCalculator(input,
							output->GetRequestedRegion(),
							gaussianOperator.GetRadius() );
		for ( fit=faceList.begin(); fit != faceList.end(); ++fit )
		{
			it = NeighborhoodIteratorType( gaussianOperator.GetRadius(),
											input, *fit );
			out = IteratorType( output, *fit );
			for (it.GoToBegin(), out.GoToBegin(); ! it.IsAtEnd(); ++it, ++out)
			{
				out.Set( innerProduct(it, gaussianOperator) );
			}
		}

		// Swap the input and output buffers
		if (d != ImageType::ImageDimension - 1)
		{
			ImageType::Pointer tmp = input;
			input = output;
			output = tmp;
		}
	}


	//
	// Writer
	//
	typedef DcmWriter< ImageType, Dimension > MyDCMWriterType;
	MyDCMWriterType writer = MyDCMWriterType(  );
	writer.SetFileName( argv[argcOutputImageFile] );
	writer.SetInput( output );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err ){
		std::cerr << "Error while writing the output image with "
			<< argv[argcOutputImageFile] << std::endl;
	}
	std::cout << "Output image wrote: " 
			<< argv[argcOutputImageFile] << std::endl;
	std::cout<<std::endl;

}//end main
