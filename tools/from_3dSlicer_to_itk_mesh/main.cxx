
/*
 * Transform a Mesh save as STL to VTK mesh that can be read with ITK
 */

#include "itkQuadEdgeMesh.h"
#include "itkSTLMeshIOFactory.h"
#include "itkSTLMeshIO.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"

// displaying results
#include "../console_tools/color.h"

int main( int argc, char* argv[]){
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
  if( argc < 2 )
	{
			error( "Missing Parameters ");
		std::cerr << "Usage: InputFileName OutputFileName" << std::endl;
		exit(EXIT_FAILURE);
	}

	blueMessage( "!!!!!!!!!! WARNING !!!!!!!!!!");
	blueMessage(std::string("Be sure that your STL mesh is save with the")
			+ std::string(" non-compress option in 3D Slicer !") );

	const unsigned int Dimension = 3;
	typedef float PixelType ;
	typedef itk :: QuadEdgeMesh < PixelType , Dimension > QEMeshType ;
	itk :: STLMeshIOFactory :: RegisterOneFactory ();


	typedef itk::MeshFileReader< QEMeshType >     ReaderType;
  typedef itk::MeshFileWriter< QEMeshType >     WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );

	writer->SetFileTypeAsASCII();

	reader->Update();
  QEMeshType * mesh = reader->GetOutput();

  //mesh->Print( std::cout );

  

  writer->SetInput( reader->GetOutput() );

  int result = EXIT_SUCCESS;

  try
  {
    writer->Update();
		finalMessage( "Output VTK mesh as: " + std::string(argv[2]) );
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << excp << std::endl;
    result = EXIT_FAILURE;
  }


	exit(result);
}



