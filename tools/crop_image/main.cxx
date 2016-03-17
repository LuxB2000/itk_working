/*
 * Crop a 3D volume based on a mesh.
 * The crop volume's sizes are based on the bounding box of the mesh.
 * INPUTS
 * ------
 * * Input image path
 * * Input Mesh defining the ImageType::Region
 * * Output image path
 * OUTPUTS
 * -------
 * * If no error crop version is written on HD
 *
 */


// my file readers
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../meshtools/meshtools.h"

#include "../console_tools/color.h"

// general
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkExtractImageFilter.h"

// mesh
#include "itkMesh.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkMeshFileReader.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 4 )
    {
    error("Missing Parameters ");
    std::cerr << "Usage:" << argv[0];
    std::cerr << " InputImagePath InputMesh outputCropImagePath"<<std::endl;
    return EXIT_FAILURE;
    }
	unsigned int argvInputImage = 1, argvInputMesh=2, argvOutputImage = 3;


	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
	typedef float						MeshPixelType;
  typedef  float          PixelType;

	// image types
  typedef itk::Image< PixelType, Dimension > ImageType;
	ImageType::Pointer image;
	ImageType::Pointer cropImage;
	typedef DcmReader<ImageType, Dimension>    MyDCMReaderType;
	typedef DcmWriter<ImageType, Dimension>    MyDCMWriterType;
	MyDCMWriterType writer = MyDCMWriterType();

	// mesh types
	typedef itk::Mesh<MeshPixelType, Dimension >   MeshType;
	typedef itk::BinaryMask3DMeshSource< ImageType, MeshType >
                                                 MeshSourceType;
	typedef itk::MeshFileReader< MeshType >				MeshReaderType;
	MeshReaderType::Pointer meshReader = MeshReaderType::New();
	//MeshSourceType::Pointer meshSource = MeshSourceType::New();


	//
	// read the mask image with my reader
	//
	MyDCMReaderType dcmReader = MyDCMReaderType(argv[argvInputImage]);
	try{
		dcmReader.Update();
		image = dcmReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		std::cerr << "Error while reading the input image with path: "<<
			argv[argvInputImage] << "." << std::endl;
		exit(EXIT_FAILURE);
	}

	//
	// Read the mesh
	//
	MeshType::Pointer mesh;
	meshReader->SetFileName( argv[argvInputMesh] );
	mesh = meshReader->GetOutput();
	try{
		meshReader->Update();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while reader the mesh with path: "<<
			argv[argvInputMesh] << "." << std::endl;
		std::cerr << err << std::endl;
	}

	//
	// Extract the region of interest
	//
	typedef MeshTools< ImageType, MeshType > MeshToolsType;
	MeshToolsType meshTools = MeshToolsType( image, mesh);
	ImageType::RegionType desiredRegion = meshTools.ComputeMeshRegion();


	//
	// crop image
	//
	typedef itk::ExtractImageFilter< ImageType, ImageType > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetExtractionRegion(desiredRegion);
  filter->SetInput(image);
#if ITK_VERSION_MAJOR >= 4
  filter->SetDirectionCollapseToIdentity(); // This is required.
#endif
  filter->Update();

	std::cout << "Output image region: " << filter->GetOutput()->GetRequestedRegion() << std::endl;


	//
	// write crop volume
	//
	writer.SetFileName( argv[argvOutputImage] );
	writer.SetInput( filter->GetOutput() );
	try{
		writer.Update();
		finalMessage( "Crop volume wrote as: "
				+ std::string(argv[argvOutputImage]) );
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while writing the crop volume as: "
			<< argv[argvOutputImage] << std::endl;
		std::cerr << err << std::endl;
	}
	
	std::cout << std::endl;
}

