/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description: Parse an input folder and set together all the data
 *								present in it in a 4 Dimensions volume
 *
 *        Version:  1.0
 *        Created:  24/11/14 13:35:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */


// my file readers
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
// displaying results
#include "../console_tools/color.h"

#include "../volume4dwriter/volume4dwriter.h"
#include "../volume4dreader/volume4dreader.h"

#include "../concentrationreader/concentrationreader.h"

// general
#include "itkImage.h"
#include "itkVectorImage.h"

// folder parser
#include <sys/types.h>
#include <dirent.h>
#include <vector>

#include "itkImageRegionIterator.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 3 )
    {
    error("Missing Parameters ");
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputFolderPath output4DImagePath removingOutliers [default:1]"<<std::endl;
    return EXIT_FAILURE;
    }
	unsigned int argvInputFolder = 1, argvOutputImage = 2;

	std::cout << "Input folder to parse: " << argv[argvInputFolder] 
		<< std::endl;

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
    return(EXIT_FAILURE);
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

	const unsigned int nmbrOfFile = fileToRead.size();

	// alphabetical ordering 
	// WARNING: assumption that alphabetical order is time ordering as well
	std::sort( fileToRead.rbegin(),
			fileToRead.rend(),
			std::greater<std::string>());

	for( int i=0; i<fileToRead.size() ; i++){
		std::cout << fileToRead.at(i) << std::endl;
	}

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  double          PixelType;

	// image types
	typedef float ArrayValType;
	typedef itk::VariableLengthVector<ArrayValType> ArrayType;
	//typedef itk::Array< ArrayValType > ArrayType;
	//typedef itk::Image< ArrayType, Dimension > ImageArrayType;
	typedef itk::VectorImage<ArrayValType,Dimension> ImageArrayType;

	typedef itk::ImageRegionIterator< ImageArrayType > IteratorArrayType;
  typedef itk::Image< PixelType, Dimension > ImageType;
	ImageType::Pointer image;
	ImageType::Pointer cropImage;
	typedef itk::ImageRegionIterator< ImageType > IteratorType;
	//typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef ConcentrationReader< ImageType, Dimension> 
																				ConcentrationReaderType;

	typedef Volume4dWriter<ImageArrayType, Dimension > My4dWriterType;
	My4dWriterType writer = My4dWriterType( );

	//ArrayType emptyArray( nmbrOfFile );
	//emptyArray.Fill( -1 );
	ArrayType emptyArray;
	emptyArray.SetSize( nmbrOfFile );
	emptyArray.Fill(0);
	ImageArrayType::Pointer volume4D = ImageArrayType::New();

	for( unsigned int f=0; f<nmbrOfFile; f++ ){
		//
		// read the mask image with my reader
		//
		//MyDCMReaderType* reader = new MyDCMReaderType();
		ConcentrationReaderType reader = ConcentrationReaderType();
		std::string path = std::string( std::string(argv[argvInputFolder]) +
				"/" + fileToRead.at(f) );
		//std::cout << "(" << f << ") " << fileToRead.at(f) << std::endl;
		reader.SetFileName( path.c_str() );
		try{
			//dcmReader->Update();
			reader.Update();
			image = reader.GetOutput();/* dcmReader->GetOutput(); */
		}
		catch( itk::ExceptionObject & err ) {
			exit(EXIT_FAILURE);
		}
	
		//
		// Initialize the 4D structure
		//
		if( f == 0 ){
			volume4D->SetOrigin( image->GetOrigin() );
			volume4D->SetRegions( image->GetRequestedRegion() );
			volume4D->SetSpacing( image->GetSpacing() );
			volume4D->SetNumberOfComponentsPerPixel( nmbrOfFile );
			try{
				volume4D->Allocate();
			}catch( itk::ExceptionObject &err ){
				std::cerr << "Error while allocated the output 4D data.\n"
					<< err
					<< std::endl;
				exit( EXIT_FAILURE );
			}
			// init the structur with empty array
			IteratorArrayType iterator( volume4D,
					volume4D->GetRequestedRegion() );
			for ( iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator ){
				iterator.Set( emptyArray );
			}
		}


		// Parse each file and set the value inside the 4D data
		IteratorType itConc( image, image->GetRequestedRegion() );
		IteratorArrayType it4D( volume4D, volume4D->GetRequestedRegion() );
		PixelType c;
		ArrayType a;
		for( itConc.GoToBegin(), it4D.GoToBegin(); !itConc.IsAtEnd(); 
				++itConc, ++it4D ){
			c = itConc.Get();
			a = it4D.Get();
			a[f] = (ArrayValType) c;
			it4D.Set( a );
			/*if( f==(nmbrOfFile-1))
				std::cout << c <<",";
			*/
		}
		std::cout << std::endl;
	}//end for each file to read

	//
	// write the 4D volume
	//
	writer.SetFileName( argv[argvOutputImage] );
	writer.SetInput( volume4D );
	try{
		writer.Update();
		finalMessage( "output 4D image in " +
				std::string( argv[argvOutputImage]));
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while the crop volume as: "
			<< argv[argvOutputImage] << std::endl;
	}


	std::cout << std::endl;
	exit(EXIT_SUCCESS);
}

