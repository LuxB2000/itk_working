/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description: Apply a deformation field to a moving image.
 *								You need to specfigy a reference image to extract the
 *								final origin and 
 *
 *        Version:  1.0
 *        Created:  18/10/14 10:52:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

/*
 */


// my file readers
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
// displaying results
#include "../console_tools/color.h"

#include "itkImageFileReader.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

// general
#include "itkImage.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 5 )
    {
    error( "Missing Parameters " );
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath ReferenceImage "
			<< "DeformationFieldPath outputImagePath"<<std::endl;
    return EXIT_FAILURE;
    }
	unsigned int argvInputImage = 1, argvReferenceImage=2,
							 argvDF=3, argvOutputImage = 4;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image, refImage;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;
	// deformation field type
	typedef itk::Vector< float, Dimension > VectorPixelType;
	typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldType;
	DisplacementFieldType::Pointer df;
	typedef itk::ImageFileReader< DisplacementFieldType >
																							DisplacementFieldReaderType;
	typedef itk::WarpImageFilter< ImageType, ImageType, DisplacementFieldType>
																											WarpImageFilterType;
	//typedef itk::LinearInterpolateImageFunction< ImageType, double>
	//																												InterpolateType;
	typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double>
																													InterpolateType;

	// read the images
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
	dcmReader2.SetFileName( argv[argvReferenceImage] );
	try{
		dcmReader2.Update();
		refImage = dcmReader2.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the reference image with: "
			<< argv[argvReferenceImage] <<std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	// read the deformation field
	DisplacementFieldReaderType::Pointer dfreader = DisplacementFieldReaderType::New();
	dfreader->SetFileName( argv[argvDF] );
	try{
		dfreader->Update();
		df = dfreader->GetOutput();
	}catch(itk::ExceptionObject err){
		std::cerr << "ERROR while reading the deformation field with: "
			<< argv[argvDF] 
			<< std::endl;
		std::cerr << err << std::endl;
	}

	// deformed the input volume
	WarpImageFilterType::Pointer warper = WarpImageFilterType::New();
	InterpolateType::Pointer interpolator = InterpolateType::New();
	warper->SetInput( image );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( refImage->GetSpacing() );
	warper->SetOutputOrigin( refImage->GetOrigin() );
	warper->SetOutputDirection( refImage->GetDirection() );
	warper->SetDisplacementField( df );
	warper->Update();

	// write the deformed volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argvOutputImage] );
	writer.SetInput( warper->GetOutput() );
	try{
		writer.Update();
		finalMessage("Write deformed image as: "
			+std::string( argv[argvOutputImage] ));
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while the volume as: "
			<< argv[argvOutputImage] << std::endl;
	}
	
	std::cout << std::endl;
	exit(EXIT_SUCCESS);
}

