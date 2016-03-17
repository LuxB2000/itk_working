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
    std::cerr << " outputImagefile zCorrection yCorrection";
		std::cerr << std::endl;
    return EXIT_FAILURE;
    }

	int argvTiffImageFile=1, argvOutputImageFile=2, 
			argvZSlice=3, argvYCoord=4;

	if( argc == 4 )
	{
		zArtefact = true;
		std::cout << "Z modif: " << atoi(argv[argvZSlice]) << std::endl;
	}else if (argc == 5 )
	{
		yArtefact = true;
	}



	//
	// main type def
	//
  const unsigned int    Dimension = 3;

  typedef  short					 OutputPixelType;

  typedef itk::Image< OutputPixelType, Dimension >  OutputImageType;

	OutputImageType::Pointer outputImage;


	//
	// Readers
	//
	typedef DcmReader<OutputImageType,Dimension> TiffDCMReaderType;
	TiffDCMReaderType tiffDCMReader =
													 TiffDCMReaderType(argv[argvTiffImageFile]);

	try{
//		std::cout<<"Read the tiff image with path: ";
//		std::cout<<argv[argvTiffImageFile]<<std::endl;
		tiffDCMReader.Update();
		outputImage = tiffDCMReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr<<"Error while reading the tiff image."<<std::endl;
		std::cerr<<err<<std::endl;
		exit(EXIT_FAILURE);
	}

	//std::cout<<outputImage<<std::endl;
	//exit(1);

	//
	// Correct Artefacts
	//
		OutputImageType::SizeType size =
                         outputImage->GetLargestPossibleRegion().GetSize();
	OutputImageType::Pointer outputImage_tmp, outputImage_tmp2;
		typedef itk::ImageRegionIterator< OutputImageType > IteratorType;
		typedef itk::ImageRegionConstIterator< OutputImageType > 
                                                         ConstIteratorType;

	if( zArtefact ){
		outputImage_tmp = OutputImageType::New();


		unsigned int zSlice = atoi( argv[argvZSlice] );

		if( (zSlice >= size[2]) || (zSlice<0) )
		{
			std::cerr<<"Error: Z slice parameter is miss formulated: "<<zSlice
				<<std::endl;
			exit(EXIT_FAILURE);
		}

		// zSlice == 0 means we don't want it
		if( zSlice != 0 ){
			outputImage_tmp->SetOrigin( outputImage->GetOrigin() );
			outputImage_tmp->SetSpacing( outputImage->GetSpacing() );
			outputImage_tmp->SetRegions( outputImage->GetLargestPossibleRegion() );
			try
			{
				outputImage_tmp->Allocate();
			}catch( itk::ExceptionObject &err )
			{
				std::cerr<<"Error while allocating the temporary image."<<std::endl;
				std::cerr<<err<<std::endl;
			}

			OutputImageType::RegionType firstRegion, secondRegion;
			OutputImageType::SizeType firstSize, secondSize;
			OutputImageType::IndexType firstStart, secondStart;

			firstStart[0] = 0;
			firstStart[1] = 0;
			firstStart[2] = zSlice;
			firstSize[0] = size[0];
			firstSize[1] = size[1];
			firstSize[2] = size[2] - zSlice; // parse zSlice to size[2]
			firstRegion.SetSize(firstSize);
			firstRegion.SetIndex(firstStart);

			secondStart[0] = 0;
			secondStart[1] = 0;
			secondStart[2] = 0;
			secondSize[0] = size[0];
			secondSize[1] = size[1];
			secondSize[2] = zSlice; // parse 0 to zSlice -1
			secondRegion.SetSize(secondSize);
			secondRegion.SetIndex(secondStart);
			std::cout << " ------------------------ " << std::endl;
			std::cout << "First Region: " << firstRegion << "Second Region: "<<secondRegion<< std::endl;

			ConstIteratorType firstConstIterator( outputImage, firstRegion ),
												secondConstIterator( outputImage, secondRegion);
			IteratorType firstIterator( outputImage_tmp, secondRegion ),
									 secondIterator( outputImage_tmp, firstRegion );

			for ( firstConstIterator.GoToBegin(), firstIterator.GoToBegin();
					!firstConstIterator.IsAtEnd(),!firstIterator.IsAtEnd();
					++firstConstIterator, ++firstIterator )
			{
				firstIterator.Set( firstConstIterator.Get() );
			}

			for ( secondConstIterator.GoToBegin(), secondIterator.GoToBegin();
					!secondConstIterator.IsAtEnd(),!secondIterator.IsAtEnd();
					++secondConstIterator, ++secondIterator )
			{
				secondIterator.Set( secondConstIterator.Get() );
			}

			outputImage = outputImage_tmp ;

		}// end if zSlice == 0
	}//end if zArtefact

	if( yArtefact ){
		unsigned int yCoord = (unsigned int) atoi(argv[argvYCoord]);
		if( yCoord != 0 ){
			outputImage_tmp = OutputImageType::New();
			outputImage_tmp->SetOrigin( outputImage->GetOrigin() );
			outputImage_tmp->SetSpacing( outputImage->GetSpacing() );
			outputImage_tmp->SetRegions( outputImage->GetLargestPossibleRegion() );
			try
			{
				outputImage_tmp->Allocate();
			}catch( itk::ExceptionObject &err )
			{
				std::cerr<<"Error while allocating the temporary image."<<std::endl;
				std::cerr<<err<<std::endl;
			}

			OutputImageType::RegionType firstRegion, secondRegion;
			OutputImageType::SizeType firstSize, secondSize;
			OutputImageType::IndexType firstStart, secondStart;


			firstStart[0] = 0;
			firstStart[1] = 0;
			firstStart[2] = 0;
			firstSize[0] = size[0];
			firstSize[1] = size[1]-yCoord;
			firstSize[2] = size[2]; // parse zSlice to size[2]
			firstRegion.SetSize(firstSize);
			firstRegion.SetIndex(firstStart);

			secondStart[0] = 0;
			secondStart[1] = yCoord;
			secondStart[2] = 0;
			secondSize[0] = size[0];
			secondSize[1] = size[1]-yCoord;
			secondSize[2] = size[2]; // parse 0 to zSlice -1
			secondRegion.SetSize(secondSize);
			secondRegion.SetIndex(secondStart);

			std::cout<< "TEST-FR " << firstRegion << std::endl;
			std::cout<< "TEST-SR " << secondRegion<< std::endl;

			ConstIteratorType firstConstIterator( outputImage, secondRegion );
			IteratorType firstIterator( outputImage_tmp, firstRegion );

			for ( firstConstIterator.GoToBegin(), firstIterator.GoToBegin();
					!firstConstIterator.IsAtEnd(),!firstIterator.IsAtEnd();
					++firstConstIterator, ++firstIterator )
			{
				firstIterator.Set( firstConstIterator.Get() );
			}

			std::cout<< "TEST-3" << std::endl;

			firstStart[0] = 0;
			firstStart[1] = size[1] - yCoord;
			firstStart[2] = 0;
			firstSize[0] = size[0];
			firstSize[1] = yCoord;
			firstSize[2] = size[2]; // parse zSlice to size[2]
			firstRegion.SetSize(firstSize);
			firstRegion.SetIndex(firstStart);

			secondStart[0] = 0;
			secondStart[1] = 0;
			secondStart[2] = 0;
			secondSize[0] = size[0];
			secondSize[1] = yCoord;
			secondSize[2] = size[2]; // parse 0 to zSlice -1
			secondRegion.SetSize(secondSize);
			secondRegion.SetIndex(secondStart);

			ConstIteratorType secondConstIterator( outputImage, secondRegion);
			IteratorType secondIterator( outputImage_tmp, firstRegion );

			for ( secondConstIterator.GoToBegin(), secondIterator.GoToBegin();
					!secondConstIterator.IsAtEnd(),!secondIterator.IsAtEnd();
					++secondConstIterator, ++secondIterator )
			{
				secondIterator.Set( secondConstIterator.Get() );
			}

			std::cout<< "TEST-BIS" << std::endl;
			outputImage = outputImage_tmp ;

		}
	}

	//
	// writer
	//
	/*
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	WriterType::Pointer writer = WriterType::New();

  writer->SetFileName( argv[argvOutputImageFile] );
	writer->SetInput( outputImage );

	try
		{
			std::cout<<"Write the image as: ";
			std::cout<<argv[argvOutputImageFile]<<std::endl;
			writer->Update();
		}
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught during the ";
		std::cerr << "image writing ! Proposed filename is ";
		std::cerr << argv[argvOutputImageFile]<<"." << std::endl; 
    std::cerr << err << std::endl; 
		}

	std::cout<<"TIF correctly transform to META and write as:";
	std::cout<<argv[argvOutputImageFile]<<std::endl;
	*/
	typedef DcmWriter< OutputImageType, Dimension > MyDCMWriterType;
	MyDCMWriterType writer = MyDCMWriterType(  );
	writer.SetFileName( argv[argvOutputImageFile] );
	writer.SetInput( outputImage );
	writer.Update();
	std::cout<<std::endl;

}//end main
