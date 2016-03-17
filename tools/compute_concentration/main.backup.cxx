/*
 * Compute the concentration of the Contrast Agent on a Post-Injection 
 * volume. We assume that the POST and PRE image are registered in the
 * same physical space.
 * We parse the PRE image, for each MESHi we find all the pixels inside the
 * MESHi. For each pixel, we find the 3D physical coordinates 
 *
 * INPUTS
 * ------
 *  * A PRE-injection image path
 *  * A POST-injection image path
 *  * A path to write 3D volume that contains concentration
 *  * A path to reach the mesh that represents the Reference signal
 *  * A list of 3D Mesh refer to the PRE image
 *  TODO: take the blood vessel signal into consideration
 *
 * OUTPUTS
 * -------
 *  * A volume of the same size of PRE. All the pixels that belong to no
 *  MESHi are set to 0. Otherwise, the concentration is computed for all the
 *  pixels belong to a MESHi
 *
 * EXAMPLE
 * -------
 *./compute preImage.mha postImage.mha outputConcentration.mha refSignal.vtk mesh1.vtk mesh2.vtk
 *
 *
 * 2014
 * Author: J Plumat, UoA
 *
 */

// my file readers
#include "../dcmreader/dcmreader.h"
#include "../concentrationwriter/concentrationwriter.h"

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


int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0]<<" --- "<<std::endl;
	std::cout << "TEST argc: " << argc << std::endl;
	for( int i=1; i<argc; i++ )
	{
		std::cout << i << " - " << argv[i] << std::endl;
	}
	// check parameters
  if( argc < 7 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " prePHANTOMImageFile preINNEREARImage"
			<< " postPHANTOMImageFile postINNEREARImage"
			<< " concentrationFilePath"
			<< " phantomMesh" <<" listOfMesh "<<std::endl;
    std::cerr << "The list must contains at least one mesh" <<std::endl;
    return EXIT_FAILURE;
    }
	const unsigned int argvPrePhantom=1, argvPreIE=2, 
				argvPostPhantom=3, argvPostIE=4,
				argvOutput=5, argvRefMesh=6,
				argvMeshListBegin = 7, argvMeshListEnd=argc;
	const unsigned int nbrOfMesh= argvMeshListEnd - argvMeshListBegin;

	std::cout << nbrOfMesh << " MESHi provided as inputs." << std::endl;

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
																	ConstIteratorType;

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
		exit(EXIT_FAILURE);
	}
	MyDCMReaderType preIEReader = MyDCMReaderType(argv[argvPreIE]);
	try{
		preIEReader.Update();
		preIE = preIEReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		exit(EXIT_FAILURE);
	}

	MyDCMReaderType postPhantomReader = MyDCMReaderType(argv[argvPostPhantom]);
	try{
		postPhantomReader.Update();
		postPhantom = postPhantomReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		exit(EXIT_FAILURE);
	}

	MyDCMReaderType postIEReader = MyDCMReaderType(argv[argvPostIE]);
	try{
		postIEReader.Update();
		postIE = postIEReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		exit(EXIT_FAILURE);
	}

	//
	// Create a volume to compute the concentration
	//
	//std::cout << preIE << std::endl;
	//std::cout << argv[argvPreIE] << std::endl;
	//std::cout << postIE << std::endl;
	concentration->SetRegions( preIE->GetLargestPossibleRegion() );
	concentration->SetSpacing( preIE->GetSpacing() );
	concentration->SetOrigin( preIE->GetOrigin() );
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
	//concentration->Fill(0.0);
	//concentration->Print( std::cout );

	// Pharmacocinetik model
	PharmaModel pharma = PharmaModel();


	//
	// Parse the reference mesh in Pre and Post and compute the correction
	// coefficients
	//
	refMeshReader = MeshReaderType::New();
	refMeshReader->SetFileName( argv[ argvRefMesh ] );
	try{
		refMeshReader->Update();
		refMesh = refMeshReader->GetOutput();
	}catch( itk::ExceptionObject &err )
	{
		std::cerr << "Error while reading the Reference Mesh as "
			<< argv[argvRefMesh] << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	typedef MRICorrections< ImageType, Dimension> MRICorrectionsType;
	MRICorrectionsType corrections = MRICorrectionsType();
	//std::cout << "test, pre phantom region: " << prePhantom->GetLargestPossibleRegion() << "origin " << prePhantom->GetOrigin() << std::endl;
	//corrections.SetPreImage( prePhantom );
	//corrections.SetPostPhantomImage( postPhantom );
	//corrections.SetPostIEImage( postIE );
	//corrections.SetReferenceMesh( refMesh );
	//corrections.Update();
	corrections.SetPostImage( postPhantom );
	corrections.SetPreImage( prePhantom );
	corrections.Update();
	MRICorrectionsType::OutputVectorType cor_vector =
		MRICorrectionsType::OutputVectorType(corrections.GetCorrectionVector());

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


	//
	// Parse each mesh
	//
	int n=0;
	for( n=0; n<nbrOfMesh; n++ )
	{
		meshReader = MeshReaderType::New();
		meshReader->SetFileName( argv[argvMeshListBegin + n] );
		std::cout << "Read the mesh: " << argv[argvMeshListBegin + n] << std::endl;
		try{
			meshReader->Update();
			mesh = meshReader->GetOutput();
		}catch( itk::ExceptionObject &err )
		{
			std::cout << "Error while reading " << n << "th mesh located at "
				<< argv[argvMeshListBegin + n ] <<std::endl;
			return EXIT_FAILURE; // TODO: more intelligent exit, ignore the mesh?
		}

		typedef MeshTools<ImageType, MeshType> MeshToolsType;
		MeshToolsType meshTools = MeshToolsType( preIE, mesh );
		ImageType::RegionType regionToParse = meshTools.ComputeMeshRegion();
		ImageType::IndexType zeros = {0,0,0};
		regionToParse.SetIndex( zeros );
		// const iterator
		//std::cout << "TEST - region:\nmesh region: " << regionToParse << "\n"
		//	"image region: " << preIE->GetLargestPossibleRegion() << " spacing: " << preIE->GetSpacing() << "origine: " << preIE->GetOrigin() <<"\n"
		//	"BB: " << mesh->GetBoundingBox()->GetBounds() << std::endl;
		ConstIteratorType preIEIterator( preIE, regionToParse );

		MeshSpatialObjectType::Pointer spatialMesh =
																MeshSpatialObjectType::New();
		spatialMesh->SetMesh( mesh );

		// Parse the PRE image
		// find all the pixels that are inside the current mesh
		// for each pixels compute its physical coordinates
		// based on the coord, find the corresponding pixels values in POST
		std::cout << ">> Go for Iterator " << std::endl;
		ImageType::IndexType preInd, postInd, concInd;
		ImageType::PointType prePoint, postPoint;
		ImageType::PixelType prePix, postPix;
		ConcentrationImageType::PixelType concPix;
		for(preIEIterator.GoToBegin(); !preIEIterator.IsAtEnd();++preIEIterator)
		{
			prePix = preIEIterator.Get();
			preInd = preIEIterator.GetIndex();
			preIE->TransformIndexToPhysicalPoint( preInd, prePoint );
			//std::cout<< " -- Index : " << preInd << " prePoint " << prePoint <<std::endl;
			if( spatialMesh->IsInside( prePoint ) )
			{
				//std::cout<<" Index: " << preInd << " prePoint: " << prePoint <<std::endl;
				postIE->TransformPhysicalPointToIndex( prePoint, postInd );
				concentration->TransformPhysicalPointToIndex( prePoint, concInd );
				postPix = postIE->GetPixel( postInd );
				concPix = (ConcentrationPixelType) pharma.ComputeConcentration( 
													(double) prePix, (double) postPix, 
													cor_vector.at(postInd[2]) );
				
				concPix = (concPix > 1E5) ? 0 : concPix;
				concPix = (concPix < 0 ) ? 0 : concPix;
				concentration->SetPixel( concInd, concPix );
			}// end if IsInside Mesh
		}// end for preIterator

	}//end for each mesh

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
    return EXIT_FAILURE;
  }

	std::cout << "Concentration volume created without any error"
		<< " and write as "<<argv[argvOutput]<<std::endl;


	//concentration->Print( std::cout );

	std::cout<<std::endl;

	exit(EXIT_SUCCESS);

}
