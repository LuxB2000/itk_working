/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *        Created:  21/10/14 18:02:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

// my file readers
#include "../dcmwriter/dcmwriter.h"
#include "../dcmreader/dcmreader.h"
#include "../volume4dreader/volume4dreader.h"
#include "../console_tools/color.h"

#include "../meshtools/meshtools.h"

// general
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageIteratorWithIndex.h"
#include "itkImageConstIteratorWithIndex.h"

// folder parser
#include <sys/types.h>
#include <dirent.h>
#include <vector>
#include <math.h>

// mesh
#include "itkMesh.h"
#include "itkMeshFileReader.h"
#include "itkMeshSpatialObject.h"
#include <iostream>

// EM
#include "itkVector.h"
#include "itkGaussianMixtureModelComponent.h"
#include "itkExpectationMaximizationMixtureModelEstimator.h"
#include "itkListSample.h"


/*
 * main
 */
int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 3 )
    {
    error( "Missing Parameters" );
    std::cerr << "Usage:" << argv[0]
			<< " referenceImage MapsFolder"
			<<std::endl;
    return EXIT_FAILURE;
    }
	unsigned int argvInputRefIm = 1, argvInputFolder = 2;

	// Constant
	std::string additional_path = "_EM-modif";
	unsigned int numberOfClasses = 2;

	std::string msg;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
	typedef  float          ArrayValType;
  typedef  float          PixelType;
  typedef  float          MeshPixelType;

	// image types
	typedef itk::VectorImage<ArrayValType,Dimension>		ImageArrayType;
  typedef itk::Image< PixelType, Dimension >					ImageType;

	ImageArrayType::RegionType region;
	ImageType::Pointer refImage, image;

	typedef itk::Vector< double, 1 >											ItkVectorType;
	typedef itk::ImageRegionIterator< ImageArrayType >		IteratorArrayType;
	typedef itk::ImageRegionIterator< ImageArrayType >		IteratorType;
	typedef itk::ImageIteratorWithIndex< ImageArrayType >
																												IteratorIndexType;
	// one sample per ROI in the reference volume
	typedef itk::Statistics::ListSample< ItkVectorType >	SampleType;
	typedef DcmReader< ImageType, Dimension >							DcmReaderType;
	typedef DcmWriter< ImageType, Dimension >							DcmWriterType;

	DcmReaderType reader, refReader;
	
	// compute statistics along each ROIs
	typedef std::vector< double >				VectorDType;
	typedef std::vector< VectorDType >	VectorVDType;
	typedef std::vector< PixelType >		RoiIDVectorType;
	
	//one sample per ROI in the reference image
	typedef std::vector< SampleType::Pointer >			SampleVectorType;
	typedef itk::Statistics::ExpectationMaximizationMixtureModelEstimator<
															SampleType > EstimatorType;
	typedef itk::Statistics::GaussianMixtureModelComponent< SampleType >
																														ComponentType;

	// contains ID of the ROI for each file in the folder
	RoiIDVectorType roiValues = RoiIDVectorType();

	//
	// Read the Reference Image
	//
	refReader = DcmReaderType();
	refReader.SetFileName( argv[argvInputRefIm] );
	try{
		refReader.Update();
		refImage = refReader.GetOutput();
		refImage->DisconnectPipeline();
	}catch(itk::ExceptionObject &err ){
		std::cerr << "Error while reading the input 4D volume with path: "
			<< argv[argvInputRefIm] << std::endl;
		std::cerr << err << std::endl;
		exit( EXIT_FAILURE );
	}

	//
	// parse folder and find all "mha" volume
	//
	struct dirent *de=NULL;
  DIR *d=NULL;
	d=opendir(argv[argvInputFolder]);
  if(d == NULL)
  {
		std::cerr << "Couldn't open input directory " <<
			argv[argvInputFolder] << std::endl;
    exit(EXIT_FAILURE);
  }
	// Loop while not NULL
	std::string fileName, ext = std::string("mha");
	std::vector<std::string> fileToRead;
  while(de = readdir(d)){
		fileName= de->d_name;
		if( fileName.substr(fileName.find_last_of(".") +1 ) == ext ){
			fileToRead.push_back( fileName );
		}
	}

	const unsigned int nbrOfFile = fileToRead.size();

	// alphabetical order
	std::sort( fileToRead.rbegin(),
			fileToRead.rend(),
			std::greater<std::string>());

	std::cout << "Input binary masks: ";
	for( int i=0; i<fileToRead.size() ; i++){
		if( i!=0 ) std::cout << ", ";
		std::cout << fileToRead.at(i);
	}
	std::cout << std::endl;

	//
	// Initial parameters of the classes `object' and `background'
	// TODO
	typedef itk::Array< double > ParametersType;
	ParametersType params(2); // mean, std and proportions
	const unsigned int nbrOfClass = 2;
	std::vector< ParametersType > initialParameters( nbrOfClass );
	// class bkg, enhancement around 0
	params[0] = 500; // estimation of mean
	params[1] = 500; // superior of std
	initialParameters[0] = params;
	// class obj, enhancement sup to 1
	params[0] = 3000; // estimation of mean
	params[1] = 500; // superior of std
	initialParameters[1] = params;

	itk::Array<double> initialProportions(nbrOfClass);
	initialProportions[0] = 0.3;
	initialProportions[1] = 0.7;


	//
	// MAIN LOOP: parse each file, find ROIs and extract mean and std
	// ---------
	bool ok = true;
	for( unsigned int f=0; f<nbrOfFile; f ++ ){
		std::cout << std::endl;
		//
		// Read each binary mask
		//
		std::cout << std::endl;
		std::string path = std::string( std::string(argv[argvInputFolder]) +
				"/" + fileToRead.at(f) );
		reader = DcmReaderType();
		reader.SetFileName( path.c_str() );
		try{
			reader.Update();
			std::cout << "Reading file: " << fileToRead.at(f) << std::endl;
			image = reader.GetOutput();
			image->DisconnectPipeline();
		}
		catch( itk::ExceptionObject & err ) {
			exit(EXIT_FAILURE);
		}

		// reset the roiValues
		roiValues = RoiIDVectorType();
		
		//
		// Find how many regions are in the atlas
		//
		ImageArrayType::IndexType ind;
		ImageType::PixelType pix;
		ImageType::PixelType roiVal;
		typedef itk::ImageRegionIteratorWithIndex< ImageType > 
																	IteratorType;
		// iterator parsing the BM
		IteratorType it( image, image->GetLargestPossibleRegion() );
		for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
			if( it.Get() != 0 ){
				// check if the vector already contains the value
				if((std::find(roiValues.begin(), roiValues.end(),it.Get()) == 
						roiValues.end()) ) {
					roiValues.push_back( it.Get() );
				}
			}
		}

		std::cout << "We find the following ROIs: ";
		for( std::vector<PixelType>::const_iterator i = roiValues.begin();
				i != roiValues.end(); ++i){
			std::cout << *i << ' ';
		}
		std::cout << std::endl;


		// initiate data structures
		SampleVectorType sample = SampleVectorType();
		for( unsigned int i=0; i<roiValues.size(); i++ ){
			SampleType::Pointer s = SampleType::New();
			s->SetMeasurementVectorSize(1); // length of meas vector
																				// ie: pixels size is 1
			sample.push_back( s );
		}

		//
		// collect samples
		//
		int c=0;
		for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
			// at site in the region
			if( it.Get() != 0 ){
				for( c=0; c<roiValues.size(); c++ ){
					// find the correct vector place
					if(roiValues.at(c)==it.Get()){
						pix = refImage->GetPixel( it.GetIndex() );
						//std::cout << pix << ", ";
						sample.at(c)->PushBack( pix );
					}
				}
			}
		}
		std::cout << std::endl;

		// EM for each sample (one sample per ROI)
		for(unsigned int r=0; r<roiValues.size(); r++){

			std::cout << ">> Roi: " << r << "\n" <<
				" sample size: " << sample.at(r)->Size() << std::endl;

				// init each class
				std::vector< ComponentType::Pointer > components;
				for( unsigned int c=0; c<nbrOfClass; c++ ){
					components.push_back( ComponentType::New() );
					(components[c])->SetSample( sample.at(r) );
					(components[c])->SetParameters( initialParameters[c] );
				} //end for l

				EstimatorType::Pointer estimator = EstimatorType::New();
				estimator->SetSample( sample.at(r) );
				estimator->SetMaximumIteration( 50 ); // TODO create variable
				estimator->SetInitialProportions( initialProportions );

				// std::cout << sample.at(r) << std::endl;

				for( unsigned int c=0; c<nbrOfClass; c++ )
				{
					estimator->AddComponent( (ComponentType::Superclass*)
						(components[c]).GetPointer() );
				}// end for c{
				try{
					estimator->Update();
				}catch( itk::ExceptionObject &err ){
					std::cerr << "Error during the EM estimation" << std::endl;
					std::cerr << err << std::endl;
					ok=false; // interrupt for the current ROIs since EM didn't work
									// and components still contains the init parameters
				}
	
				// extract parameters
				if( ok ){ // the EM ends well
					std::vector< ParametersType > classParameters( nbrOfClass );
					itk::Array<double> proportions(nbrOfClass);
					for ( unsigned int c = 0; c < nbrOfClass; c++ )
					{
						// mean and std
						params[0] = (components[c])->GetFullParameters()[0]; 
						params[1] = sqrt( (components[c])->GetFullParameters()[1] );
						classParameters[c] = params;
						proportions[c] = estimator->GetProportions()[c];

						std::cout << "mean: " << params[0] << 
							" std: " << params[1] << 
							" proportion: " << proportions[c] <<
							std::endl;
					}

					// modify the ROI in the map
					// hyp: the class with the highest value ( should be the second)
					// is the obj, every pixel with value > or < to mean+/-std should
					// be considered as outside of the class
					long N =0;
					for( it.GoToBegin(); !it.IsAtEnd(); ++it ){
						if( it.Get() == roiValues.at(r) ){
							pix = refImage->GetPixel( it.GetIndex() );
							if( 
								pix > abs( classParameters[1][0]+1*classParameters[1][1] ) 
								|| 
								pix < abs( classParameters[1][0]-1*classParameters[1][1] ) 
								){
									// remove the index from the map
									N ++;
									it.Set( 0 );
							}
						}
					}

					std::cout << "We Remove " << N << " points." << std::endl;

				} //end if ok

				ok = true;
				
		}//end for r
	
		// save the new maps
		//std::cout << fileToRead.at(f) << std::endl;
		std::size_t found = fileToRead.at(f).find_last_of(".");
		std::string newName ;
		if (found!=std::string::npos){
			newName = fileToRead.at(f).substr(0,found)
				+ additional_path + fileToRead.at(f).substr(found);
			std::string newPath =
				std::string( std::string(argv[argvInputFolder]) +
				"/" + newName );
			std::cout << newName << std::endl;
		}
		DcmWriterType dcmwriter = DcmWriterType();
		path = std::string( std::string(argv[argvInputFolder]) +
				"/" + newName );
		dcmwriter.SetFileName( path.c_str() );
		dcmwriter.SetInput( image );
		try{
			dcmwriter.Update();
			std::cout << "Modified maps written as " << newName << std::endl;
		}catch( itk::ExceptionObject &err ){
			std::cerr << "Error while writing the modified maps with "
				<< path << std::endl;
			exit( EXIT_FAILURE );
		}

	} // end for each file

	finalMessage( "ROIs map modified.");

	std::cout << std::endl;
}

