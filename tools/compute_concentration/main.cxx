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

// pharma models
#include "pharma.h"

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

// deal with writing results
std::ofstream fileStream;
void closeFile(){
	fileStream.close();
}

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0]<<" --- "<<std::endl;
	std::cout << "TEST argc: " << argc << std::endl;
	for( int i=1; i<argc; i++ )
	{
		std::cout << i << " - " << argv[i] << std::endl;
	}
	// check parameters
  if( argc < 6 )
  {
		error( "Missing Parameters ");
		std::cerr << "Usage: " << argv[0];
		std::cerr << " prePHANTOMImageFile preINNEREARImage"
			<< " postPHANTOMImageFile postINNEREARImage"
			<< " concentrationFilePath\n" 
			<< "[PhantomValueThresh=1000] "
			<< "[FileNameToWriteCorrCoef=/tmp/results.txt]" << std::endl;
		return EXIT_FAILURE;
  }
	const unsigned int argvPrePhantom=1, argvPreIE=2, 
				argvPostPhantom=3, argvPostIE=4,
				argvOutput=5,argvThresh=6, argcFileName=7;

	double thresh = 1000;
	if( argc > argvThresh ){
		thresh = atof(argv[argvThresh]);
	}
	std::cout << "Threshold value: " << thresh << std::endl;

	bool writeCoef = false;
	std::string fileName = std::string("/tmp/results.txt");
	if( argc > argcFileName ){
		fileName = std::string( argv[argcFileName] );
	}
	// write the information at the end of 
	fileStream.open( fileName.c_str(), std::ofstream::app );
	if( !fileStream.is_open() ){
		std::cerr << "Impossible to write corr coefficient in file: "
			<< fileName << std::endl;
		return EXIT_FAILURE;
	}else{
		writeCoef = true;
		std::cout << "Writing the correction coefficient in "
			<< fileName << std::endl;
	}

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

	ImageType::Pointer prePhantom, preIE;
	ImageType::Pointer postPhantom, postIE;
	ConcentrationImageType::Pointer concentration = 
				ConcentrationImageType::New();

	MeshType::Pointer mesh, refMesh;
	MeshReaderType::Pointer meshReader, refMeshReader;

	//
	// read the pre and post image with my reader
	//
	typedef DcmReader<ImageType,Dimension> MyDCMReaderType;
	MyDCMReaderType preReader = MyDCMReaderType(argv[argvPrePhantom]);
	try{
		preReader.Update();
		prePhantom = preReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		closeFile();
		std::cerr << "Error while reading pre phantom: " <<
			argv[argvPrePhantom] << std::endl;
		exit(EXIT_FAILURE);
	}
	MyDCMReaderType preIEReader = MyDCMReaderType(argv[argvPreIE]);
	try{
		preIEReader.Update();
		preIE = preIEReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		closeFile();
		std::cerr << "Error while reading pre IE: " <<
			argv[argvPreIE] << std::endl;
		exit(EXIT_FAILURE);
	}

	MyDCMReaderType postPhantomReader = MyDCMReaderType(argv[argvPostPhantom]);
	try{
		postPhantomReader.Update();
		postPhantom = postPhantomReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		closeFile();
		std::cerr << "Error while reading post phantom: " <<
			argv[argvPostPhantom] << std::endl;
		exit(EXIT_FAILURE);
	}

	MyDCMReaderType postIEReader = MyDCMReaderType(argv[argvPostIE]);
	try{
		postIEReader.Update();
		postIE = postIEReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		closeFile();
		std::cerr << "Error while reading post IE: " <<
			argv[argvPostIE] << std::endl;
		exit(EXIT_FAILURE);
	}

	//
	// Create a volume to compute the concentration
	//
	concentration->SetRegions( postIE->GetLargestPossibleRegion() );
	concentration->SetSpacing( postIE->GetSpacing() );
	concentration->SetOrigin( postIE->GetOrigin() );
	try 
  { 
    concentration->Allocate();
  } 
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught during the concentration ";
		std::cerr <<"volume memory allocation!"<< std::endl; 
    std::cerr << err << std::endl; 
		closeFile();
    return EXIT_FAILURE;
  }

	std::cout << " -- TEST 0 -- " << std::endl;

	// Pharmacocinetik model
	PharmaModel pharma = PharmaModel();
	double pixelVol = 1.0;
	for( unsigned int i=0; i<3; i++){
		pixelVol = pixelVol * postIE->GetSpacing()[i];
	}

	std::cout << " -- TEST 1 -- " << std::endl;


	typedef MRICorrections< ImageType, Dimension> MRICorrectionsType;
	MRICorrectionsType corrections = MRICorrectionsType();
	corrections.SetPostImage( postPhantom );
	corrections.SetPreImage( prePhantom );
	corrections.SetThresholdValue( thresh );
	corrections.Update();
	std::cout << " -- TEST 2 -- " << std::endl;
	MRICorrectionsType::OutputVectorType cor_vector =
		MRICorrectionsType::OutputVectorType(corrections.GetCorrectionVector());
	MRICorrectionsType::OutputVectorType pre_vector =
		MRICorrectionsType::OutputVectorType(corrections.GetPreCoeff());
	MRICorrectionsType::OutputVectorType noize_pre_mean_vector =
		MRICorrectionsType::OutputVectorType(corrections.GetPreNoizeMean());
	MRICorrectionsType::OutputVectorType noize_pre_std_vector =
		MRICorrectionsType::OutputVectorType(corrections.GetPreNoizeStd());
	MRICorrectionsType::OutputVectorType noize_post_mean =
		MRICorrectionsType::OutputVectorType(corrections.GetPostNoizeMean());
	MRICorrectionsType::OutputVectorType noize_post_std =
		MRICorrectionsType::OutputVectorType(corrections.GetPostNoizeStd());

	std::cout<< "Correction coefficients (one per slice): [";
	int l=cor_vector.size();
	for(int i=0;i<l;i++){
		if(i!=0){
			std::cout<<", ";
		}
		std::cout << "(" << i <<")"
						<< cor_vector.at(i);
	}
	std::cout<<"]"<<std::endl;

	// write coefficient in file
	fileStream << "PreVolumePath," << argv[argvPreIE] << "\n";
	fileStream << "PostVolumePath," << argv[argvPostIE] << "\n";
	fileStream << "PostRatio," << corrections.GetPostRatio() << "\n";
	fileStream << "PreRatio," << corrections.GetPreRatio() << "\n";
	fileStream << "Threshold," << thresh << "\n";
	fileStream << "Slice"
		<< ",Corrective coefficients"
		<<",Pre coef only"
		<<",noize pre mean,noize pre std"
		<<",noize post mean,noize post std"
		<<"\n";

	for(int i=0;i<l;i++){
		fileStream << i <<",";
		fileStream << cor_vector.at(i) << ",";
		fileStream << pre_vector.at(i) << "," ;
		fileStream << noize_pre_mean_vector.at(i)<<",";
		fileStream << noize_pre_std_vector.at(i) << "," ;
		fileStream << noize_post_mean.at(i)<<","<<noize_post_std.at(i);
		fileStream << "\n" ;
	}

	fileStream << "\n --------------------------------------------- \n" ;
	
	IteratorWithIndexType itPre( preIE,
														preIE->GetLargestPossibleRegion().GetSize() );
	ConcentrationImageType::PixelType concPix;
	ImageType::PixelType prePix, postPix, phantomPix;
	ImageType::IndexType preInd, postInd, concInd, phantomInd;
	ImageType::PointType prePoint;

	double corr_coef = 1.0;
	bool is_inside = false;

	for( itPre.GoToBegin(); !itPre.IsAtEnd(); ++itPre ){
		preInd = itPre.GetIndex();
		prePix = itPre.Get();
		preIE->TransformIndexToPhysicalPoint( preInd, prePoint );
		
		is_inside = postIE->TransformPhysicalPointToIndex( prePoint, postInd );
		postPix = postIE->GetPixel( postInd );

		if( preInd[2] != postInd[2] ){
			std::cout << " >>\t " << preInd << " --- " << postInd <<std::endl;
		}

		prePhantom->TransformPhysicalPointToIndex( prePoint, phantomInd );
//		if( is_inside ){
			corr_coef = cor_vector.at(preInd[2]); //phantomInd[2]);
//		}else{
//			corr_coef = 1;
//			std::cout << " --- strange " << std::endl;
//		}

		

		concPix = (ConcentrationPixelType)
							pharma.ComputeConcentration(
											(double) prePix, (double) postPix,
											corr_coef, pixelVol
							);
		concentration->TransformPhysicalPointToIndex( prePoint, concInd );

		//std::cout<< prePix<< " " << postPix << " " << cor_vector.at(preInd[2])
		//	<< " " << concPix << std::endl;

		//std::cout << "PrePix=" << prePix << " PostPix=" << postPix << std::endl;

		double before = concPix ;
		concPix = (concPix > 1E5) ? 0 : concPix;
		concPix = (concPix < 0 ) ? 0 : concPix;

		/*
		if( concPix == 0 ){
			std::cout << concInd << " " << before << " pre:" << prePix << " post:" << postPix << " corr:" << corr_coef
				<< " isInside? " << is_inside << " preInd:" << preInd << " postInd:" << postInd << std::endl;
		}
		*/

		concentration->SetPixel( concInd, concPix );
	}

	//
	// Write the concentration volume
	//
	//caster->SetInput( concentration );
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
		closeFile();
    exit( EXIT_FAILURE );
  }

	finalMessage(
			"Concentration volume created without any error and write as "
		 + std::string(argv[argvOutput]) );


	//concentration->Print( std::cout );

	std::cout<<std::endl;

	closeFile();
	exit(EXIT_SUCCESS);

}
