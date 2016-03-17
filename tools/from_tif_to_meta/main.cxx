/*
 * Transform a tiff image and save it as a META image
 *
 * INPUTS
 * * a path to reach the TIFF image,
 * * a path to save the META image.
 * * an artifact along the Z axis that needs to be corrected [optinal,
 * default: none]. This parameter is a single value. It designates a slice
 * number zSlice. After correction, all the slices from zSlice to size[3] 
 * are moved as first slices of the volumes and zSlice is now the first of 
 * them. Notice: the order is not changed. Set it to 0 if you don't want to 
 * change.
 * * an artifact along the Y axis that needs to be corrected [optinal,
 * default: none] the coordinate designate the space between size[2] and
 * Ycoord
 *
 * OUTPUTS
 * * An image is save on the HD.
 *
 * 2014
 * Author: J Plumat
 *
 */

// my dcm reader
#include "./../dcmreader/dcmreader.h"
#include "./../dcmwriter/dcmwriter.h"

// general
#include "itkImage.h"

// iterators
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"


/*
 * MAIN
 */
int main( int argc, char *argv[] )
{
	std::cout<<" -- "<<argv[0]<<" -- "<<std::endl;
	bool zArtefact = false, yArtefact = false;
	// check parameters
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " TiffImageFile";
    std::cerr << " outputImagefile [zCorrection, default=0]";
		std::cerr << std::endl;
    return EXIT_FAILURE;
    }

	int argvTiffImageFile=1, argvOutputImageFile=2, 
			argvZSlice=3;

	//
	// main type def
	//
  const unsigned int    Dimension = 3;
  typedef float					 OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension >  ImageType;
	typedef itk::ImageRegionIterator< ImageType > IteratorType;

	ImageType::Pointer outputImage;


	//
	// Readers
	//
	typedef DcmReader<ImageType,Dimension> TiffDCMReaderType;
	TiffDCMReaderType reader =
													 TiffDCMReaderType(argv[argvTiffImageFile]);

	try{
		reader.Update();
		outputImage = reader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr<<"Error while reading the input image."<<std::endl;
		std::cerr<<err<<std::endl;
		exit(EXIT_FAILURE);
	}


	//
	// Allocate a tmp volume
	//
	ImageType::Pointer tmpImage = ImageType::New();
	tmpImage->SetOrigin( outputImage->GetOrigin() );
	tmpImage->SetSpacing( outputImage->GetSpacing() );
	tmpImage->SetRegions( outputImage->GetRequestedRegion() );
	try{
		tmpImage->Allocate();
	}catch( itk::ExceptionObject &err ){
		std::cerr << "Error while allocating a temporary image."
			<< err << std::endl;
	}

	//
	// Correct artefact
	//
	unsigned int zSlice = 0;
	if( argc==4 )
	{
		zSlice = atoi( argv[argvZSlice] );
	}
	std::cout << "Correction slice factor: " << zSlice << std::endl;

	ImageType::SizeType sz1, sz2, sz3, sz4;
	ImageType::IndexType ind1, ind2;
	ImageType::RegionType from0toZslice,fromZsclietoEnd;

	sz1 = outputImage->GetRequestedRegion().GetSize();
	sz1[2] = zSlice;
	ind1.Fill(0);
	from0toZslice.SetIndex( ind1 );
	from0toZslice.SetSize( sz1 );

	sz2 = outputImage->GetRequestedRegion().GetSize();
	ind2.Fill(0);
	ind2[2] = sz2[2] - zSlice;
	sz2[2] = zSlice;
	fromZsclietoEnd.SetIndex( ind2 );
	fromZsclietoEnd.SetSize( sz2 );

	//std::cout << "from0toZslice: " << from0toZslice 
	//	<< "\n fromZsclietoEnd: " << fromZsclietoEnd << std::endl;

	IteratorType it(outputImage, from0toZslice );
	IteratorType it2(tmpImage, fromZsclietoEnd );
	for( it.GoToBegin(), it2.GoToBegin(); !it.IsAtEnd(), !it2.IsAtEnd();
			++it, ++it2 ){
		it2.Set( it.Get() );
	}

	// --- 
	sz1 = outputImage->GetRequestedRegion().GetSize();
	sz1[2] = sz1[2]-zSlice;
	ind1.Fill(0);
	ind1[2] = zSlice;
	from0toZslice.SetIndex( ind1 );
	from0toZslice.SetSize( sz1 );

	sz2 = outputImage->GetRequestedRegion().GetSize();
	ind2.Fill(0);
	sz2[2] = sz2[2] - zSlice;
	fromZsclietoEnd.SetIndex( ind2 );
	fromZsclietoEnd.SetSize( sz2 );

	//std::cout << "from0toZslice: " << from0toZslice 
	//	<< "\n fromZsclietoEnd: " << fromZsclietoEnd << std::endl;

	IteratorType it3(outputImage, from0toZslice );
	IteratorType it4(tmpImage, fromZsclietoEnd );
	for( it3.GoToBegin(), it4.GoToBegin(); !it3.IsAtEnd(), !it4.IsAtEnd();
			++it3, ++it4 ){
		it4.Set( it3.Get() );
	}

	//
	// Writer
	//
	typedef DcmWriter< ImageType, Dimension > MyDCMWriterType;
	MyDCMWriterType writer = MyDCMWriterType(  );
	writer.SetFileName( argv[argvOutputImageFile] );
	writer.SetInput( tmpImage );
	writer.Update();
	std::cout<<std::endl;

}//end main
