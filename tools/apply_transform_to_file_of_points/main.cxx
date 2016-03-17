/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description: read a input file of points and an input transform. It applies the 
 *    transform to the points and write the result in an output file.
 *    WARNING: the transform is suppose to be 3DVesor.
 *    TODO improve by making the transform more flexible
 *    see http://itk-users.7.n7.nabble.com/ITK-users-applying-transform-to-point-td35198.html
 *    Inputs
 *			- file to read index
 *			- transform
 *			- output file
 *
 *
 *        Version:  1.0
 *        Created:  23/04/15 10:22:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */



// my file readers
#include "../file_tools/readIndexFile.h"
#include "../file_tools/writeIndexFile.h"
// displaying results
#include "../console_tools/color.h"

// mesh reader and writer
#include "itkQuadEdgeMesh.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"
#include "itkTransformMeshFilter.h"

// general
#include <itkTransformFileReader.h>
#include "itkSimilarity3DTransform.h"
#include "itkAffineTransform.h"

#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImage.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
	const unsigned int argvInputFile = 1, argvInputTransform = 2,
							 argvOutputFile = 3;
  if( argc < 3 )
    {
    error( "Missing Parameters ");
    std::cerr << "Usage:" << argv[0] << " ";
    std::cerr << "InputFile InputTransform"
			<< " outputFile";
		std::cerr << std::endl;
    return EXIT_FAILURE;
    }

	//
	// main type def
	//
	bool verbose = true;
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;
	 typedef itk::NearestNeighborInterpolateImageFunction<
												ImageType, double > InterpolatorType;

	// image types

	// mesh type
	

	// transform
	typedef itk::TransformFileReaderTemplate<double>::TransformListType
		TransformListType;
	TransformListType* transformList;
  //typedef itk::AffineTransform< double, Dimension> TransformType;
  //typedef itk::Similarity3DTransform< double> TransformType;
  typedef itk::VersorRigid3DTransform< double> TransformType;
  typedef itk::ResampleImageFilter< 
                            ImageType, 
                            ImageType >       ResampleFilterType;
	ResampleFilterType::Pointer resampler;

	// read the input file
	typedef ReadIndexFile<ImageType::PointType, Dimension> ReadIndexType;
	ReadIndexType indexReader = ReadIndexType();
	indexReader.SetFileName( argv[argvInputFile] );
	indexReader.Update();
	ReadIndexType::VectorIndexType pointList = ReadIndexType::VectorIndexType( indexReader.GetIndexList() );
	unsigned int nbrOfPoints = pointList.size();
	std::cout << nbrOfPoints << " found in the file." << std::endl;


	//
	// Read the input transform
	//
#if (ITK_VERSION_MAJOR==4 && ITK_VERSION_MINOR>=5) || ITK_VERSION_MAJOR>4
	itk::TransformFileReaderTemplate<double>::Pointer transformReader = 
		itk::TransformFileReaderTemplate<double>::New();
#else
  itk::TransformFileReader::Pointer transformReader =
		itk::TransformFileReader::New();
#endif
	transformReader->SetFileName( argv[argvInputTransform] );
	try{
		transformReader->Update();
		transformList = transformReader->GetTransformList();
	}catch( itk::ExceptionObject &err ){
		std::cerr << "Error while reading the input transfom: " 
			<<argv[argvInputTransform]<< "\n" << err << std::endl;
		exit( EXIT_FAILURE );
	}

	// For each transform
	TransformListType::iterator it;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	for( it=transformList->begin(); it!=transformList->end(); ++it )
	{
		TransformType::Pointer affineTransform = 
									static_cast<TransformType*> ( (*it).GetPointer() );
		resampler = ResampleFilterType::New();

		//std::cout << affineTransform <<std::endl;

		// deform each point
		for( unsigned int i=0; i<nbrOfPoints; i++ ){
			std::cout << pointList.at(i) << " -> " ;
			pointList.at(i) = affineTransform->TransformPoint( pointList.at(i) );
			std::cout << pointList.at(i) << std::endl;
		}

	}

	// write the deformed index
	typedef WriteIndexFile<ImageType::PointType, Dimension> WriteIndexType;
	WriteIndexType indexWriter = WriteIndexType();
	indexWriter.SetFileName( argv[argvOutputFile] );
	indexWriter.SetIndexList( pointList );
	indexWriter.Update();

	std::cout << std::endl;
	exit(EXIT_SUCCESS);
}

