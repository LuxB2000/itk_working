/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description: Create a rigid transform based on rotation (in degrees) and translation.
 *							The transformation will be saved in file.
 *
 *        Version:  1.0
 *        Created:  20/03/16 09:05:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

// my libraries
#include "../dcmreader/dcmreader.h"
#include "../dcmwriter/dcmwriter.h"
#include "../console_tools/color.h"

#include "itkImage.h"
#include "itkTransformFileWriter.h"

#include "itkRigid3DTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "itkAffineTransform.h"

int main( int argc, char *argv[] )
{
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 8 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
    std::cerr << " rotLR rotIS rotPA transX transY transZ outputTransformPath (rationAtOrigin)"<<std::endl;
    return EXIT_FAILURE;
  }
	unsigned int argvRLR = 1, argvRIS = 2, argvRPA = 3,
							 argvTX = 4,  argvTY = 5,  argvTZ = 6,
							 argvOutputPath = 7,
							 argcrotationAtCenter = 8,
							 argcIsCenter = 9;

	//
	// main type def
	//
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef  itk::Image<PixelType, Dimension> ImageType;

	// image types
	ImageType::Pointer image;
	//typedef itk::Rigid3DTransform< float > TransformType;
	typedef itk::AffineTransform< double > TransformType;
	//typedef itk::VersorRigid3DTransform< double > TransformType;

  itk::TransformFileWriterTemplate<double>::Pointer transformWriter =
			itk::TransformFileWriterTemplate<double>::New();

	float tx = atof( argv[argvTX] ),
	ty = atof( argv[argvTY] ),
	tz = atof( argv[argvTZ] ),
	// angles in radians
	rlr = atof( argv[argvRLR] ) * atan(1.0f)/45.0,
	ris = atof( argv[argvRIS] ) * atan(1.0f)/45.0,
	rla = atof( argv[argvRPA] ) * atan(1.0f)/45.0;

  TransformType::Pointer outputTransform = TransformType::New();
	
	// set rotations
	outputTransform->Rotate( 1, 2, rlr );
	outputTransform->Rotate( 0, 2, ris );
	outputTransform->Rotate( 0, 1, rla );

	//TransformType::ParametersType params[outputTransform->GetNumberOfParameters()];
	//TransformType::OffsetType::ValueType offTransflationV[3] = {tx,ty,tz};
	//TransformType::OffsetType offTranslation = offTransflationV;

//	outputTransform->SetOffset( offTranslation );

	//TODO: 
	// the center of rotation is by default the origin
	// outputTransform->SetCenter();
	//outputTransform->Rotate(rlr, ris, rla);
	
	/*
	TransformType::ParametersType param;
	param.fill(0);

	param[10] = tx;
	param[11] = ty;
	param[12] = tz;

	param[3] = rlr;
	param[4] = ris;
	param[5] = rla;

	//outputTransform->SetParameters( param );

TransformType::OutputVectorType transl; 
	transl[0] = tx;
	transl[1] = ty;
	transl[2] = tz;
	//outputTransform->SetTranslation( transl );
	//
	//
*/


	transformWriter->SetFileName( argv[argvOutputPath] );
	transformWriter->SetInput( outputTransform );
	try{
			transformWriter->Update();
	}catch(itk::ExceptionObject &err ){
		std::cerr << "Error while writing the transform as: ";
		std::cerr << argv[argvOutputPath] << "\n";
		std::cerr << err << std::endl;
	}

	std::string msg = "Transformation wrote with T=(" + SSTR(tx)
		+ ", " + SSTR(ty) + ", " + SSTR(tz) + ") and R=(" + SSTR(rlr)
			+ ", " + SSTR(ris) + " , " + SSTR(rla) + ") saved at ";
		msg += argv[argvOutputPath];
	finalMessage( msg );
	std::cout << std::endl;
}

