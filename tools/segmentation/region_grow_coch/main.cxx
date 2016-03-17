/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  Use the region growing to segment the cochlea. This pipeline use an
 *    edge detector (sigmoid of the gradient magnitude) to isolate the cochlea. Based on a 
 *    set of seed (specified in input) the region growing finds part of the cochlea and 
 *    outputs a binary map. This is improved with a closing.
 *    Inputs:
 *			- input and output paths, 
 *			- sigma for smoothing (0.05) for GP cochlea segmentation
 *			- alpha and beta for sigmoid parameterization (-300 and 2000)
 *			- upper and lower threshold for binarization (0 and 0.3)
 *			- radius for closing ball (1)
 *			- path to read the file with coordinates of the seeds
 *
 *
 *    THIS PIPELINE REQUIRES A USER VERIFICATION
 *
 *        Version:  1.0
 *        Created:  29/04/15 15:19:33
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */


//#include "itkImageFileReader.h"
//#include "itkImageFileWriter.h"
#include "../../dcmreader/dcmreader.h"
#include "../../dcmwriter/dcmwriter.h"
#include "../../console_tools/color.h"
#include "../../file_tools/readIndexFile.h"
#include "../../file_tools/writeIndexFile.h"
//#include "../../console_tools/debug.h"

#include <itkCastImageFilter.h>
#include "itkGrayscaleFunctionDilateImageFilter.h"

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"

#include "itkConnectedThresholdImageFilter.h"

#include <itkClosingByReconstructionImageFilter.h>
#include <itkBinaryClosingByReconstructionImageFilter.h>
#include "itkBinaryBallStructuringElement.h"


int main( int argc, char *argv[] )
{
	std::cout << " ---- " << argv[0] << " ---- " << std::endl;

	///////////////////////
	// Deal with inputs
	///////////////////////
  typedef   float           InternalPixelType;
  const     unsigned int    Dimension = 3;
  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
	InternalImageType::IndexType seedInd;

	const int argcInputFileName = 1, argcOutputFileName = 2,
						argcSigma = 3,
						argcAlpha = 4, argcBeta = 5,
						argcUppererThresh = 6, argcLowerThresh = 7,
						argcRaduis = 8,
						argcSeedpath = 9;

	std::string msg = "";

	if( argc <= argcSeedpath ){
		error( "Missing parameters");
		std::cerr << "[InputFileName] [OutputFileName] [Sigma] [Alpha] [Beta] "
			<< "[upperThers] [lowerThres] [radius] [fileToReadSeed]"
			<< std::endl;
			exit( EXIT_FAILURE );
	}

	// extract seeds
	unsigned int nbrOfSeeds = 0;
	/*
	if( ( (argc-argcBeginSeed) % 3 ) != 0 ){
		error( "You should enter a list of 3 coordinates for each seeds." );
		exit( EXIT_FAILURE );
	}

	nbrOfSeeds = (argc-argcBeginSeed) / 3;
	msg = SSTR( nbrOfSeeds ) + std::string( " seed(s) enter: " );
	*/

	unsigned int seedC=0, s=0, i=0;
	//for( unsigned int s=0; s++; s<nbrOfSeeds ){
	/*
	while( s<nbrOfSeeds  ){
		i=0;
		while( i<Dimension ){
			seedC = atoi( argv[ argcBeginSeed + (s)*(Dimension) + i ] );
			seedInd[i] = seedC;
	//		std::cout << (s*Dimension) << " + " << (i) << " = " << argcBeginSeed + (s*Dimension)+(i) << std::endl;
			i ++;
		}
		seedList.push_back( seedInd );
		s++;
	}
	*/
	typedef ReadIndexFile<InternalImageType::PointType, Dimension> ReadIndexType;
	ReadIndexType pointReader = ReadIndexType();
	pointReader.SetFileName( argv[argcSeedpath] );
	pointReader.Update();
	ReadIndexType::VectorIndexType seedList = ReadIndexType::VectorIndexType( pointReader.GetIndexList() );


	//

	nbrOfSeeds = seedList.size();

	//for( unsigned int s=0; s++; s<nbrOfSeeds ){
	s = 0;
	msg = SSTR( nbrOfSeeds ) + std::string( " seed(s) enter: " );
	while( s<nbrOfSeeds  ){
		msg += std::string("[") +  SSTR( seedList.at(s)[0] ) + std::string(", ") +
					SSTR( seedList.at(s)[1] ) + std::string(", ") +
					SSTR( seedList.at(s)[2] ) + std::string("] ") ;
		s++;
	}

	blueMessage( msg );
	///////////////////////
	// Main typedef
	///////////////////////
	typedef DcmReader<InternalImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<InternalImageType, Dimension>    MyDCMWriterType;

	typedef itk::ImageFileWriter<InternalImageType> DebugWriterType;

  MyDCMReaderType reader = MyDCMReaderType();
  MyDCMWriterType writer = MyDCMWriterType();
  reader.SetFileName( argv[argcInputFileName] );
	writer.SetFileName( argv[argcOutputFileName] );

  typedef   itk::CurvatureAnisotropicDiffusionImageFilter<
                               InternalImageType,
                               InternalImageType >  SmoothingFilterType;
  SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();

  typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter<
                               InternalImageType,
                               InternalImageType >  GradientFilterType;
  typedef   itk::SigmoidImageFilter<
                               InternalImageType,
                               InternalImageType >  SigmoidFilterType;
  GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();
  SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();

	typedef itk::ConnectedThresholdImageFilter< InternalImageType,
																				InternalImageType > ConnectedFilterType;
	ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();

	typedef itk::BinaryBallStructuringElement<
							InternalPixelType,
							Dimension > StructuringElementType;
	StructuringElementType structuringElement;

	typedef itk::GrayscaleFunctionDilateImageFilter<
															InternalImageType, InternalImageType,
															StructuringElementType > ClosingFilterType;
	ClosingFilterType::Pointer closingFilter = ClosingFilterType::New();


	reader.Update();
	InternalImageType::Pointer vol = reader.GetOutput();
	/////////////////////////
	// parameters
	/////////////////////////
  smoothing->SetTimeStep( 0.003 );
  smoothing->SetNumberOfIterations(  5 );
  smoothing->SetConductanceParameter( 9.0 );

  sigmoid->SetOutputMinimum(  0.0  );
  sigmoid->SetOutputMaximum(  1.0  );

  const double sigma = atof( argv[argcSigma] );
  gradientMagnitude->SetSigma(  sigma  );

  const double alpha =  atof( argv[argcAlpha] );
  const double beta  =  atof( argv[argcBeta] );
  sigmoid->SetAlpha( alpha );
  sigmoid->SetBeta(  beta  );

	float lowerThreshold = atof(argv[argcLowerThresh]);
	float upperThreshold = atof(argv[argcUppererThresh]);
	connectedThreshold->SetLower( lowerThreshold );
	connectedThreshold->SetUpper( upperThreshold );
	connectedThreshold->SetReplaceValue( 1 );
	std::cout << "ConnectedThreshold with upper=" <<
					upperThreshold << " and lower=" << lowerThreshold << std::endl;
	s = 0;
	while( s<nbrOfSeeds ){
		//seedInd = seedList.at(s);
		vol->TransformPhysicalPointToIndex( seedList.at(s), seedInd );
		connectedThreshold->SetSeed( seedInd );
		s++;
	}



	int radius = atoi( argv[argcRaduis] );
	StructuringElementType::SizeType kernelSize;
	for( unsigned int d=0; d<Dimension; d++ ){
			kernelSize[d] = radius;
	}
	structuringElement.SetRadius( kernelSize );
	structuringElement.CreateStructuringElement();
	closingFilter->SetKernel( structuringElement );
	//closingFilter->FullyConnectedOn();
	//closingFilter->SetForegroundValue( 1 );
	//structuringElement.Print(std::cout);
	std::cout << "Closing with a ball, radius=" << radius << std::endl;

	/////////////////////////
	// Pipeline
	/////////////////////////
  smoothing->SetInput( reader.GetOutput() );
  gradientMagnitude->SetInput( smoothing->GetOutput() );
  sigmoid->SetInput( gradientMagnitude->GetOutput() );
	connectedThreshold->SetInput( sigmoid->GetOutput() );
	closingFilter->SetInput( connectedThreshold->GetOutput() );
	writer.SetInput( closingFilter->GetOutput() );

	DebugWriterType::Pointer writer1 = DebugWriterType::New();
	writer1->SetInput( sigmoid->GetOutput() );
	writer1->SetFileName("/tmp/test1.mha");

	try{
		writer1->Update();
		writer.Update();
	}catch(itk::ExceptionObject &err ){
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}

	msg = std::string("Segmented choclea wrote as ") + SSTR(argv[argcOutputFileName]);
	finalMessage( msg );

	std::cout << std::endl;
	exit( EXIT_SUCCESS );
}
