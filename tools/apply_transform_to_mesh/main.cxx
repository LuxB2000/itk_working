/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description: deform a VTK mesh image with rigid or affine input
 *								transforms.
 *								WARNING: we assume that all the inputs transforms are
 *								affine transforms.
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
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
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
  if( argc < 4 )
    {
    error( "Missing Parameters ");
    std::cerr << "Usage:" << argv[0] << " ";
    std::cerr << "InputMeshPath InputTransform"
			<< " deformMeshPath";
		std::cerr << std::endl;
    return EXIT_FAILURE;
    }
	const unsigned int argvInputMesh = 1, argvInputTransform = 2,
							 argvOutputMesh = 3;

	//
	// main type def
	//
	bool verbose = true;
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;
	 typedef itk::NearestNeighborInterpolateImageFunction<
												ImageType, double > InterpolatorType;
	typedef itk :: QuadEdgeMesh < PixelType , Dimension > QEMeshType ;
	QEMeshType::Pointer mesh;

	// image types

	// mesh type
	typedef itk::MeshFileReader< QEMeshType >     ReaderType;
  typedef itk::MeshFileWriter< QEMeshType >     WriterType;
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
	

	// transform
	typedef itk::TransformFileReaderTemplate<double>::TransformListType
		TransformListType;
	TransformListType* transformList;
  typedef itk::AffineTransform< double, Dimension> TransformType;
  //typedef itk::Similarity3DTransform< double> TransformType;
  typedef itk::ResampleImageFilter< 
                            ImageType, 
                            ImageType >       ResampleFilterType;
	ResampleFilterType::Pointer resampler;
	typedef itk::TransformMeshFilter< QEMeshType, QEMeshType,
																									TransformType > TransformMeshFilterType;

	// read the input mesh
	reader->SetFileName( argv[argvInputMesh] );
	try{
		reader->Update();
		mesh = reader->GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the input mesh with: "
			<< argv[argvInputMesh] <<std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

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
	TransformMeshFilterType::Pointer transformFilter;
	TransformListType::iterator it;
	InterpolatorType::Pointer interpolator;
	QEMeshType::Pointer deformed = mesh;
	for( it=transformList->begin(); it!=transformList->end(); ++it )
	{
		TransformType::Pointer affineTransform = 
									static_cast<TransformType*> ( (*it).GetPointer() );

		if(verbose){
			std::cout << affineTransform << std::endl;
		}

		transformFilter = TransformMeshFilterType::New();
		transformFilter->SetTransform( affineTransform );
		transformFilter->SetInput( deformed );
		transformFilter->Update();
		deformed = transformFilter->GetOutput();
	}

	// write the deformed mesh
	writer->SetFileName( argv[argvOutputMesh] );
	writer->SetInput( transformFilter->GetOutput() );
	try{
		writer->Update();
		finalMessage( "Write the deformed mesh as: "
			+ std::string(argv[argvOutputMesh]) );
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while the deformed volume as: "
			<< argv[argvOutputMesh] << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}
	
	std::cout << std::endl;
	exit(EXIT_SUCCESS);
}

