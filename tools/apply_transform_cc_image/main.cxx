/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description: deform an input image with rigid or affine input
 *								transforms.
 *								WARNING: we assume that all the inputs transforms are
 *								affine transforms.
 *
 *        Version:  1.0
 *        Created:  20/10/14 09:05:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */



// my file readers
#include "../concentrationreader/concentrationreader.h"
#include "../concentrationwriter/concentrationwriter.h"
// displaying results
#include "../console_tools/color.h"

// general
#include <itkTransformFileReader.h>
#include "itkAffineTransform.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImage.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 5 )
    {
    error( "Missing Parameters ");
    std::cerr << "Usage:" << argv[0] << " ";
    std::cerr << "InputImagePath InputDeformation"
			<< " ReferenceImage deformImagePath";
		std::cerr << std::endl;
    return EXIT_FAILURE;
    }
	const unsigned int argvInputImage = 1, argvInputTransform = 2,
							 argvInputReferenceImage = 3, argvOutputImage = 4;

	//
	// main type def
	//
	bool verbose = false;
	int defaultOutputPixel = 0;
  const    unsigned int   Dimension = 3;
  typedef  double          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;
	 typedef itk::NearestNeighborInterpolateImageFunction<
												ImageType, double > InterpolatorType;

	// image types
	ImageType::Pointer image, refImage;
	typedef ConcentrationReader<ImageType, Dimension>    MyDCMReaderType;
	typedef ConcentrationWriter<ImageType, Dimension>    MyDCMWriterType;

	// transform
	typedef itk::TransformFileReaderTemplate<double>::TransformListType
		TransformListType;
	TransformListType* transformList;
  typedef itk::AffineTransform< double, Dimension> TransformType;
  typedef itk::ResampleImageFilter< 
                            ImageType, 
                            ImageType >       ResampleFilterType;
	ResampleFilterType::Pointer resampler;

	// read the image
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argvInputImage] );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the input image with: "
			<< argv[argvInputImage] <<std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	MyDCMReaderType dcmReader2 = MyDCMReaderType();
	dcmReader2.SetFileName( argv[argvInputReferenceImage] );
	try{
		dcmReader2.Update();
		refImage = dcmReader2.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the reference image with: "
			<< argv[argvInputReferenceImage] <<std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	//
	// Read the transform
	//
#if (ITK_VERSION_MAJOR==4 && ITK_VERSION_MINOR>=5) || ITK_VERSION_MAJOR>4
	itk::TransformFileReaderTemplate<double>::Pointer transformReader = 
		itk::TransformFileReaderTemplate<double>::New();
#else
  itk::TransformFileReader::Pointer transformReader =
		itk::TransformFileReader::New();
#endif
	transformReader->SetFileName( argv[argvInputTransform] );
	try{
		transformReader->Update();
		transformList = transformReader->GetTransformList();
	}catch( itk::ExceptionObject &err ){
		std::cerr << "Error while reading the input transfom: " 
			<<argv[argvInputTransform]<< "\n" << err << std::endl;
		exit( EXIT_FAILURE );
	}


	// For each transform
	TransformListType::iterator it;
	InterpolatorType::Pointer interpolator;
	ImageType::Pointer deformed = image;
	ImageType::SizeType sz = /* refImage*/ image->GetLargestPossibleRegion().GetSize() ;
	ImageType::PointType origin = refImage->GetOrigin();
	for( it=transformList->begin(); it!=transformList->end(); ++it )
	{
		interpolator = InterpolatorType::New();
		resampler = ResampleFilterType::New();
		TransformType::Pointer affineTransform = 
									static_cast<TransformType*> ( (*it).GetPointer() );
		if(verbose){
			std::cout << affineTransform << std::endl;
		}

		if(verbose){
			std::cout << origin << std::endl;
		}

		resampler->SetTransform( affineTransform );
		resampler->SetInterpolator( interpolator );
		resampler->SetInput( deformed );
		resampler->SetSize( sz );
		resampler->SetOutputOrigin(  origin );
		resampler->SetOutputSpacing( refImage->GetSpacing() );
		resampler->SetOutputDirection( refImage->GetDirection() );
		resampler->SetDefaultPixelValue( defaultOutputPixel );
		resampler->Update();
		deformed = resampler->GetOutput();
	}

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argvOutputImage] );
	writer.SetInput( resampler->GetOutput() );
	try{
		writer.Update();
		finalMessage( "Write the deformed image as: "
			+ std::string(argv[argvOutputImage]) );
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while the deformed volume as: "
			<< argv[argvOutputImage] << std::endl;
		std::cerr << err << std::endl;
	}
	
	std::cout << std::endl;
	exit(EXIT_SUCCESS);
}

