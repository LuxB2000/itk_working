/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  filter an image by preserving the edges
 *
 *        Version:  1.0
 *        Created:  09/02/15 10:30:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * The filtering is performed with a gradient magintude anistropic diffusion
 * filter.
 *
 * The maximum allowable time step for this filter is 1/2^N, where N is the 
 * dimensionality of the image. For 2D images any value below 0.250 is
 * stable, and for 3D images, any value below 0.125 is stable.
 * The conductance parameter governs the sensitivity of the conductance
 * equation. 
 * =====================================================================================
 */


// my libraries
#include "../../dcmreader/dcmreader.h"
#include "../../dcmwriter/dcmwriter.h"
#include "../../console_tools/color.h"

// general
#include "itkImage.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 2 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath outputImagePath"
			<< " [NbrOfIteration, default=5] [conductance, default=1]"
			<< " [TimeStep, default=1/2^Dimension]"
			<<std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argcInputImage = 1, argcOutputImage = 2,
		argcNbrOfIter = 3, argcConductance = 4, argcTimeStep = 5;

  const    unsigned int   Dimension = 3;

	int nbrOfIterations=5;
	float conductance = 1, timeStep=1/((float)pow(2,Dimension));

	if( argc > argcNbrOfIter ){
		nbrOfIterations = atoi( argv[argcNbrOfIter]);
	}
	if( argc > argcConductance ){
		conductance = atof( argv[argcConductance]);
	}
	if( argc > argcTimeStep ){
		timeStep = atof( argv[argcTimeStep] );
	}

	std::cout << "Filtering with parameters: nbrOfiter=" << nbrOfIterations 
		<< " conductance=" << conductance << " timeStep=" << timeStep
		<< std::endl;


	//
	// main type def
	//
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;

	// read the image
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argcInputImage] );
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();// dcmReader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	// Filtering
	typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType,
																	ImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetInput( image );
	filter->SetNumberOfIterations( nbrOfIterations );
	filter->SetTimeStep( timeStep );
	filter->SetConductanceParameter( conductance );

	// write the volume
	MyDCMWriterType writer = MyDCMWriterType( );
	writer.SetFileName( argv[argcOutputImage] );
	writer.SetInput( filter->GetOutput() );
	try{
		writer.Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the volume as: "
			<< argv[argcOutputImage] << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}
	
	std::string msg = "Output volume wrote with ";
	msg += argv[argcOutputImage];
	finalMessage( msg );
	std::cout << std::endl;
}

