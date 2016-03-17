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
#include "../concentrationreader/concentrationreader.h"
#include "../concentrationwriter/concentrationwriter.h"
// displaying results
#include "../console_tools/color.h"

// general
#include "itkImage.h"
#include "itkCastImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageFileWriter.h"

// deal with writing results
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
		std::cerr << " preImageFile postImageFile"
			<< " concentrationFilePath\n" 
			<< std::endl;
		return EXIT_FAILURE;
  }
	const unsigned int argvPreImage=1, argvPostImage=2, 
				argvOutput=3;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  double	        ConcentrationPixelType;
  typedef  double					OutputPixelType;

  typedef itk::Image< ConcentrationPixelType, Dimension > 
																						ConcentrationImageType;
  typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;

	typedef itk::ImageRegionIteratorWithIndex< ConcentrationImageType > 
																	ConcentrationIteratorType;

	typedef itk::CastImageFilter< ConcentrationImageType,
																OutputImageType> CasterType;
	CasterType::Pointer caster = CasterType::New();

	ConcentrationImageType::Pointer concentration, preNorm, postNorm;
	ConcentrationImageType::SizeType sz;

	//
	// read the pre and post image with my reader
	//
	typedef ConcentrationReader< ConcentrationImageType, Dimension> ConcentrationReaderType;
	ConcentrationReaderType readerPre = ConcentrationReaderType();
	ConcentrationReaderType readerPost = ConcentrationReaderType();

	readerPre.SetFileName( argv[argvPreImage] );
	try{
		readerPre.Update();
		preNorm = readerPre.GetOutput();
	}catch( itk::ExceptionObject &err ){
		std::cerr << "Error while reading " << argv[argvPreImage] << std::endl;
		std::cerr << err << std::endl;
		return( EXIT_FAILURE );
	}
	
	sz = preNorm->GetLargestPossibleRegion().GetSize();
	
	readerPost.SetFileName( argv[argvPostImage] );
	try{
		readerPost.Update();
		postNorm = readerPost.GetOutput();
	}catch( itk::ExceptionObject &err ){
		std::cerr << "Error while reading " << argv[argvPostImage] << std::endl;
		std::cerr << err << std::endl;
		return( EXIT_FAILURE );
	}

	/*
	for( unsigned int d=0; d<Dimension; d++ ){
		int s = postNorm->GetLargestPossibleRegion().GetSize()[d] ;
		if( s < sz[d] ){
			sz[d] = s;
		}
	}
	*/

	std::cout << "Size of image: " << sz << std::endl;
	

	//
	// Create a volume to compute the concentration
	//
	concentration = ConcentrationImageType::New();
	concentration->SetRegions( preNorm->GetLargestPossibleRegion() );
	concentration->SetSpacing( preNorm->GetSpacing() );
	concentration->SetOrigin( preNorm->GetOrigin() );
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

	//
	// Compute enhancement
	//
	ConcentrationImageType::IndexType pre_ind, post_ind;
	ConcentrationImageType::PixelType pre_pix, post_pix, conc_pix = 1;
	ConcentrationImageType::PointType pre_pt, post_pt;
	bool isInside;

	ConcentrationIteratorType it(preNorm, sz);
	for( it.GoToBegin(); !it.IsAtEnd(); ++it){
		pre_ind = it.GetIndex();
		pre_pix = it.Get();

		preNorm->TransformIndexToPhysicalPoint( pre_ind, pre_pt );
		isInside = postNorm->TransformPhysicalPointToIndex( pre_pt, post_ind );

		conc_pix = 1;

		if( isInside ){
			post_pix = postNorm->GetPixel( post_ind );
			conc_pix =  (ConcentrationPixelType) (post_pix/pre_pix);
			if( (conc_pix > 1E3) || (conc_pix < 0) ){
				conc_pix = 1;
			}
		}
		concentration->SetPixel( pre_ind, conc_pix );
	}//end for

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
    exit( EXIT_FAILURE );
  }

	finalMessage(
			"Concentration volume created without any error and write as "
		 + std::string(argv[argvOutput]) );


	//concentration->Print( std::cout );

	std::cout<<std::endl;

	exit(EXIT_SUCCESS);

}
