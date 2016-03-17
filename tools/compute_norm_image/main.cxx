/*
 * Compute the concentration of the Contrast Agent on a Post-Injection 
 * volume. We assume that the POST and PRE image are registered in the
 * same physical space.
 * We parse the PRE image, for each MESHi we find all the pixels inside the
 * MESHi. For each pixel, we find the 3D physical coordinates 
 *
 * INPUTS
 * ------
 *  * PRE-injection image paths: phantom and region of interest
 *  * POST-injection image paths: phantom and region of interest
 *  * ouput concentration volume path
 *  * [optional] threshold value, default: 1000
 *  * [optional] file to write numerical results, default /tmp/results.txt
 *				numerical results:
 *					Ph_t, Ph_0, noise (mean and std of bkg)
 *					Phantom (Ph) is defined as > thresold
 *					background (bkg) is defined as < threshold
 *  TODO: take the blood vessel signal into consideration
 *
 * OUTPUTS
 * -------
 *  coefficient in file  
 *
 * EXAMPLE
 * -------
 *
 *
 * 2014
 * Author: J Plumat, UoA
 *
 */

// my file readers
#include "../dcmreader/dcmreader.h"
#include "../concentrationwriter/concentrationwriter.h"
// displaying results
#include "../console_tools/color.h"

// dealing with MRI images
#include "MRIcorrections.h"
#include "../meshtools/meshtools.h"

// general
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageFileWriter.h"

// mesh
#include "itkMesh.h"
#include "itkMeshFileReader.h"
#include "itkMeshSpatialObject.h"

// writing results
#include <iostream>

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0]<<" --- "<<std::endl;
	std::cout << "TEST argc: " << argc << std::endl;
	for( int i=1; i<argc; i++ )
	{
		std::cout << i << " - " << argv[i] << std::endl;
	}
	// check parameters
  if( argc < 4 )
  {
		error( "Missing Parameters ");
		std::cerr << "Usage: " << argv[0];
		std::cerr
			<< " PhantomImageFile Image"
			<< " NormalisedImage\n" 
			<< "[PhantomValueThresh=1000] "
			<< std::endl;
		return EXIT_FAILURE;
  }
	const unsigned int argvPhantom=1, argvImage=2,
				argvOutput=3,argvThresh=4;

	double thresh = 1000;
	if( argc > argvThresh ){
		thresh = atof(argv[argvThresh]);
	}
	std::cout << "Threshold value: " << thresh << std::endl;

	//std::cout << nbrOfMesh << " MESHi provided as inputs." << std::endl;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  double          PixelType;
  typedef  double	        ConcentrationPixelType;
  typedef  double					OutputPixelType;
  typedef  float          MeshPixelType;

  typedef itk::Image< PixelType, Dimension >        ImageType;
  typedef itk::Image< ConcentrationPixelType, Dimension > 
																						ConcentrationImageType;
  typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;

	typedef itk::ImageRegionIterator< ConcentrationImageType > 
																	ConcentrationIteratorType;
	typedef itk::ImageRegionConstIteratorWithIndex< ImageType > 
																	ConstIteratorWithIndexType;
	typedef itk::ImageRegionIteratorWithIndex< ImageType > 
																	IteratorWithIndexType;

	typedef itk::CastImageFilter< ConcentrationImageType,
																OutputImageType> CasterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	CasterType::Pointer caster = CasterType::New();
	//WriterType::Pointer writer = WriterType::New();

	typedef itk::Mesh<MeshPixelType, Dimension >   MeshType;
	typedef itk::MeshFileReader< MeshType >				 MeshReaderType;
	typedef itk::MeshSpatialObject< MeshType >		 MeshSpatialObjectType;

	ImageType::Pointer phantom, image;

	//
	// read the pre and post image with my reader
	//
	typedef DcmReader<ImageType,Dimension> MyDCMReaderType;
	MyDCMReaderType phantomReader = MyDCMReaderType(argv[argvPhantom]);
	try{
		phantomReader.Update();
		phantom = phantomReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading phantom: " <<
			argv[argvPhantom] << std::endl;
		exit(EXIT_FAILURE);
	}
	MyDCMReaderType imageReader = MyDCMReaderType(argv[argvImage]);
	try{
		imageReader.Update();
		image = imageReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading image: " <<
			argv[argvImage] << std::endl;
		exit(EXIT_FAILURE);
	}

	//
	// Create a volume to compute the concentration
	//
	ConcentrationImageType::Pointer concentration = ConcentrationImageType::New();
	concentration->SetRegions( image->GetLargestPossibleRegion() );
	concentration->SetSpacing( image->GetSpacing() );
	concentration->SetOrigin(  image->GetOrigin() );
	try 
  { 
    concentration->Allocate();
  } 
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught during the concentration ";
		std::cerr <<"volume memory allocation!"<< std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
  }

	ImageType::SizeType sz = image->GetLargestPossibleRegion().GetSize();
	ImageType::SizeType szPhant = phantom->GetLargestPossibleRegion().GetSize();
	ImageType::IndexType ind, ind2; ind[0]=0; ind[1]=0; ind[2]=0;
	ind2[2]=0; ind2[1]=0; ind2[0]=0;
	ImageType::PixelType pix;
	ConcentrationImageType::PixelType conc;

	int minSz2=sz[2];
	if( minSz2 > szPhant[2] ) minSz2 = szPhant[2];

	long mean =0;
	long N=0;

	//
	// For each slice, compute the mean of the phantom (each pixel > Threshold)
	// Along the slice, each pixel of the image will be divided by the mean of the phantom
	//
	std::cout << "Mean along the slizes: ";
	for( ind[2]=0; ind[2]<minSz2; ind[2]++ ){

		// compute the mean along the phantom
		ind2[0]=0; ind2[1]=0; ind2[2]=ind[2];
		mean = 0;
		N =0;
		for( ind2[0]=0; ind2[0]<szPhant[0]; ind2[0]++ ){
		for( ind2[1]=0; ind2[1]<szPhant[1]; ind2[1]++ ){
			pix = phantom->GetPixel( ind2 );
			if( pix>=thresh ){
				mean += pix;
				N += 1;
			}
		}
		}

		double meanD = ((double) mean) / ((double) N);
		if( N == 0 ){
			meanD = thresh;
		}
		std::cout << "(" << ind[2] <<")" << meanD <<" ";

		// normalized the image
		for( ind[1]=0; ind[1]<sz[1]; ind[1]++ ){
		for( ind[0]=0; ind[0]<sz[0]; ind[0]++ ){
			pix = image->GetPixel( ind );
			conc = (ConcentrationImageType::PixelType) (((double) pix)/meanD);
			concentration->SetPixel( ind, conc );
		}
		}

	}
	std::cout<< std::endl;

	//
	// Write the concentration volume
	//
	typedef ConcentrationWriter< ConcentrationImageType, Dimension >
																							ConcentrationWriterType;
	ConcentrationWriterType writer = ConcentrationWriterType();
	writer.SetInput( concentration );
	writer.SetFileName( argv[argvOutput] );
	//writer->SetInput( caster->GetOutput() );
	//writer->SetFileName( argv[argvOutput] );
	try{
		writer.Update();
	}catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught during the output image writing !";
		std::cerr << std::endl; 
    std::cerr << err << std::endl; 
    exit( EXIT_FAILURE );
  }

	finalMessage(
			"Normalised volume created without any error and write as "
		 + std::string(argv[argvOutput]) );


	//concentration->Print( std::cout );

	std::cout<<std::endl;

	exit(EXIT_SUCCESS);

}
