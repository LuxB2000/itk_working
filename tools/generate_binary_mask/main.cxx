/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  Generate a 3D binary mask based on a input image and 
 *									two 3d coordinates. The two coordinates represent the 
 *									bottom left corner (blc) and the top right corner (trc).
 *									Then we can assume that blc[i] < trc[i] for all i=[0,2].
 *									The output binary mask has the sizes, origin, spacing,
 *									etc. than the input image.
 *
 *        Version:  1.0
 *        Created:  13/11/14 12:05:14
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

// general
#include "itkImage.h"

	int
main ( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] << " --- " <<std::endl;
	// check parameters
  if( argc < 9 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " InputImagePath outputBinaryMask " <<
			"blc[0] blc[1] blc[2] " <<
			"trc[0] trc[1] trc[2]"<<std::endl;
		std::cout << std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argvInputImage = 1, argvOutputMask = 2, 
							 argvblc0=3, argvblc1=4, argvblc2=5,
							 argvtrc0=6, argvtrc1=7, argvtrc2=8;


	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;

  typedef itk::Image< PixelType, Dimension >        ImageType;
	typedef ImageType::IndexType											IndexType;

	typedef itk::ImageRegionIterator< ImageType > IteratorType;

	typedef DcmReader<ImageType, Dimension> DcmReaderType;
	typedef DcmWriter<ImageType, Dimension> DcmWriterType;

	// corners
	IndexType blc, trc;
	blc[0] = atoi( argv[argvblc0] );
	blc[1] = atoi( argv[argvblc1] );
	blc[2] = atoi( argv[argvblc2] );

	trc[0] = atoi( argv[argvtrc0] );
	trc[1] = atoi( argv[argvtrc1] );
	trc[2] = atoi( argv[argvtrc2] );


	std::cout << "Input Button Left Corner " << blc << " Top Right Corner " <<
		trc << std::endl;

	// images
	ImageType::Pointer originImage;
	ImageType::Pointer maskImage = ImageType::New();

	// output binary mask pixel values
	const PixelType defaultOutputPixel = 0;
	const PixelType positiveOutputPixel = 1;

	//
	// Read the input
	//
	DcmReaderType reader = DcmReaderType();
	reader.SetFileName( argv[argvInputImage] );
	try{
		reader.Update();
		originImage = reader.GetOutput();
	}catch(itk::ExceptionObject &err ){
		std::cerr << "Error while reading the input image with path: " <<
			argv[argvInputImage] << std::endl;
		return( EXIT_FAILURE );
	}

	// ouputs some information
	ImageType::SizeType size = originImage->GetLargestPossibleRegion().GetSize();
	std::cout<<"The ouput image will have the following parameters:\n";
	std::cout<<"\t    size: "<<size<<"\n";
	std::cout<<"\t  origin: "<<originImage->GetOrigin()<<"\n";
	std::cout<<"\t spacing: "<<originImage->GetSpacing()<<"\n";

	//
	// copy information and Allocated
	//
	maskImage->SetRegions( originImage->GetLargestPossibleRegion() );
	maskImage->SetSpacing( originImage->GetSpacing());
	maskImage->SetOrigin( originImage->GetOrigin() );
	try 
  { 
    maskImage->Allocate();
  } 
  catch( itk::ExceptionObject & err ) 
  { 
    std::cerr << "ExceptionObject caught during the mask image";
		std::cerr <<" allocation!"<< std::endl; 
    std::cerr << err << std::endl; 
    return EXIT_FAILURE;
  }


	//
	// Assign the mask
	// TODO: give the possibility to keep the original image instead of
	// defaultOutputPixel
	IteratorType iterator( maskImage, maskImage->GetRequestedRegion() );
	for ( iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator )
	{
		iterator.Set( defaultOutputPixel );
	}

	//
	// assign the 1 pixel
	//
	IndexType ind; ind.Fill(0);
	for( ind[2]=blc[2]; ind[2]<=trc[2]; ind[2]++ )
		for( ind[1]=blc[1]; ind[1]<=trc[1]; ind[1]++ )
			for( ind[0]=blc[0]; ind[0]<=trc[0]; ind[0]++ )
			{
				if( ind[0]<0 || ind[1]<0 || ind[2]<0 || ind[0]>=size[0]
						|| ind[1]>=size[1] || ind[2]>=size[2] )
				{
					std::cerr<<"Warning: index out of bounds: "<<ind;
					std::cerr<<", we just ignore it."<<std::endl;
				}else{
					maskImage->SetPixel(ind, positiveOutputPixel);
				}
			}

	//
	// Write Binary Mask
	//
	DcmWriterType writer = DcmWriterType();
	writer.SetFileName( argv[argvOutputMask] );
	writer.SetInput( maskImage );
	try{
		writer.Update();
	}catch(itk::ExceptionObject &err){
		std::cerr << "Error while writing the ouput binary mask with path "<<
			argv[argvOutputMask] << std::endl;
	}
	std::cout << std::endl;
	return EXIT_SUCCESS;
}	/* ----------  end of function main  ---------- */

