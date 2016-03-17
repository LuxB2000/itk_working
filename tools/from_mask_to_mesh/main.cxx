/*
 * Create 3D Mesh based on a 3D Binary Mask and write the mesh in HD
 * INPUTS
 * ------
 * * Binary mask path
 * * Output mesh path
 * * the value of foreground pixel in binary mask (optional, default: 1)
 * * the value of the background pixel in binary mask (optional, default: 0)
 * OUTPUTS
 * -------
 * * If no error a mesh is written on HD
 *
 */


// my file readers
#include "../dcmreader/dcmreader.h"

// general
#include "itkImage.h"
#include "itkImageFileReader.h"

// mesh
#include "itkBinaryMask3DMeshSource.h"
#include "itkMesh.h"
#include "itkMeshFileWriter.h"

#include <cfloat>
// displaying results
#include "../console_tools/color.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 2 )
    {
    error( "Missing Parameters" );
    std::cerr << "Usage:" << argv[0];
    std::cerr << " binaryMaskFile outputMeshPath"<<std::endl;
		std::cerr << " foregroundValue (optional)";
		std::cerr << " backgroundValue (optional)";
		std::cerr << std::endl;
    return EXIT_FAILURE;
    }
	unsigned int argvMaskInput = 1, argvMeshOutput = 2,
							 argvForeground = 3, argvBackground = 4;


	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          MeshPixelType;
  typedef  float					MaskPixelType;

	// image types
  typedef itk::Image< MaskPixelType, Dimension > MaskImageType;
	MaskImageType::Pointer mask;
	typedef DcmReader<MaskImageType, Dimension>    MyDCMReaderType;

	// mesh types
	typedef itk::Mesh<MeshPixelType, Dimension >   MeshType;
	typedef itk::BinaryMask3DMeshSource< MaskImageType, MeshType >
                                                 MeshSourceType;
	typedef itk::MeshFileWriter< MeshType >        MeshWriterType;
	MeshWriterType::Pointer meshWriter = MeshWriterType::New();
	MeshSourceType::Pointer meshSource = MeshSourceType::New();



	MaskPixelType foreground = 1, background = 0;
	if( argc == argvBackground +1 ){
		background = static_cast<MaskPixelType>( atof( argv[argvBackground] ) );
	}
	if( argc >= argvForeground + 1 ){
		foreground = static_cast<MaskPixelType>( atof( argv[argvForeground] ) );
	}
	std::cout<<"Parse the binary mask with [foreground/background] = [";
	std::cout<<foreground<<"/"<<background<<"]"<<std::endl;
	meshSource->SetObjectValue( foreground );

	//
	// read the mask image with my reader
	//
	MyDCMReaderType dcmReader = MyDCMReaderType();
	dcmReader.SetFileName( argv[argvMaskInput] );
	try{
		dcmReader.Update();
		mask = dcmReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
		exit(EXIT_FAILURE);
	}

	//std::cout<<mask<<std::endl;
	meshSource->SetInput( mask );
	try
	{
    meshSource->Update();
  }
  catch( itk::ExceptionObject & exp )
  {
    std::cerr << "Exception thrown during Update() " << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
  }
	std::cout << "Nodes = " << meshSource->GetNumberOfNodes() << std::endl;
  std::cout << "Cells = " << meshSource->GetNumberOfCells() << std::endl;

	meshWriter->SetInput( meshSource->GetOutput() );
	meshWriter->SetFileName( argv[argvMeshOutput] );

	try
	{
		meshWriter->Update();
		finalMessage("Output wrote as: " + std::string(argv[argvMeshOutput]));
	}catch( itk::ExceptionObject &err)
	{
		std::cerr << "Error while writing the mesh with the path: ";
		std::cerr << argv[argvMeshOutput]<<std::endl;
		std::cerr << err <<std::endl;
		exit(EXIT_FAILURE);
	}

	std::cout<<"Mesh created without any error and write as ";
	std::cout<<argv[argvMeshOutput]<<std::endl;
	std::cout<<std::endl;

	if( 0 )
	{
		typedef MaskImageType::IndexType IndexType;
		typedef MeshType::PointType PointType;
		typedef MeshType::CellType CellType;
		typedef MeshType::PointsContainer::ConstIterator PointsIterator;
		PointsIterator pointsIterator = meshSource->GetOutput()->GetPoints()->Begin();
		PointsIterator pointsEnd = meshSource->GetOutput()->GetPoints()->End();
		int c=0;
		IndexType ind;
		while( c<10 && pointsIterator!= pointsEnd )
		{
			PointType p = pointsIterator.Value(); // access the point
			 mask->TransformPhysicalPointToIndex( p, ind );
			
			std::cout << "point: "<< p << " index: ";
			std::cout << ind;
			std::cout << std::endl;

			++pointsIterator;
			c++;
		}
	}


	std::cout << std::endl;
	return EXIT_SUCCESS;
}
