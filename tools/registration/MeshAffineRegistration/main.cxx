/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  Register two meshes with affine transform
 *
 *        Version:  1.0
 *        Created:  23/04/15 12:44:31
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

// my libraries
#include "../../dcmreader/dcmreader.h"
#include "../../dcmwriter/dcmwriter.h"
#include "../../console_tools/debug.h"

// mesh reader and writer
#include "itkQuadEdgeMesh.h"
#include "itkMeshFileReader.h"
#include "itkMeshFileWriter.h"
#include "itkTransformMeshFilter.h"
#include "itkTransformFileWriter.h"

// registration
#include "itkPointSetToPointSetRegistrationMethod.h"
#include "itkEuclideanDistancePointMetric.h"
#include "itkLevenbergMarquardtOptimizer.h"

#include "itkTranslationTransform.h"
#include "itkAffineTransform.h"
#include "itkEuler3DTransform.h"

#include "itkTransformMeshFilter.h"

int main( int argc, char *argv[] )
{
	// INTRO
	// -----
	std::cout<<" --- "<<argv[0] <<" --- " <<std::endl;
	// check parameters
  if( argc < 3 )
  {
    error("Missing Parameters") ;
    std::cerr << "Usage:" << argv[0];
		// TODO: more inputs
    std::cerr << " FixedMeshPath MovingMeshPath outputMeshPath"
			<< "[directTransform] [inverseTransform]"
			<< std::endl;

			exit( EXIT_FAILURE );
	}
	// ARGC
	// ------
	const unsigned int argcFixedMesh = 1, argcMovingMesh = 2, argcOutputMesh = 3,
		argcDirectTransform = 4, argcInverseTransform = 5;

	std::string msg;
	bool saveTransform = false;
	if( argc >=5 ){
		/*
		msg = std::string("Save the direct and inverse transform with\n") 
			+ argv[argcDirectTransform] + "\n" 
			+ argv[argcInverseTransform];
		debug(msg);
		*/
		saveTransform = true;
	}

	// main type def
	// --------------
	// gloabl type
  const    unsigned int   Dimension = 3;
  typedef  float          PixelType;
	typedef itk::QuadEdgeMesh < PixelType , Dimension > QEMeshType ;
	QEMeshType::Pointer fixed, moving;

	// mesh type
	typedef itk::MeshFileReader< QEMeshType >     ReaderType;
  typedef itk::MeshFileWriter< QEMeshType >     WriterType;
  WriterType::Pointer writer = WriterType::New();
  ReaderType::Pointer fReader = ReaderType::New(); // fixed reader
  ReaderType::Pointer mReader = ReaderType::New(); // moving reader

	// Load inputs
	// ------------
	fReader->SetFileName( argv[argcFixedMesh] );
	try{
		fReader->Update();
		fixed = fReader->GetOutput();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while reading the fixed mesh with "
			<< argv[argcFixedMesh] << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	mReader->SetFileName( argv[argcMovingMesh] );
	try{
		mReader->Update();
		moving = mReader->GetOutput();
	}catch( itk::ExceptionObject &err){
		std::cerr << "Error while reading the moving mesh with "
			<< argv[argcFixedMesh] << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	// configure registration
	// --------------------------
	typedef itk::PointSetToPointSetRegistrationMethod<
                                            QEMeshType,
                                            QEMeshType > RegistrationType;

	typedef itk::EuclideanDistancePointMetric<
                                    QEMeshType,
                                    QEMeshType> MetricType;

	typedef itk::AffineTransform< double, Dimension >      TransformType;
	//typedef itk::TranslationTransform< double, Dimension >      TransformType;	
  // Optimizer Type
	typedef itk::LevenbergMarquardtOptimizer OptimizerType;

	MetricType::Pointer         metric        = MetricType::New();
	TransformType::Pointer      transform     = TransformType::New();
	RegistrationType::Pointer   registration  = RegistrationType::New();
	OptimizerType::Pointer      optimizer     = OptimizerType::New();
		optimizer->SetUseCostFunctionGradient(false);
		OptimizerType::ScalesType scales( transform->GetNumberOfParameters() );

	const double translationScale = 1;   // dynamic range of translations
	const double rotationScale    = 1;   // dynamic range of rotations

	scales[0] = 1.0 / rotationScale;
	scales[1] = 1.0 / rotationScale;
	scales[2] = 1.0 / rotationScale;
	/*
	scales[3] = 1.0 / translationScale;
	scales[4] = 1.0 / translationScale;
	scales[5] = 1.0 / translationScale;
	*/
	
	unsigned long   numberOfIterations =  1000; // TODO: make it dynamic
	double          gradientTolerance  =  1e-4;   // convergence criterion
	double          valueTolerance     =  1e-4;   // convergence criterion
	double          epsilonFunction    =  1e-5;   // convergence criterion

	optimizer->SetScales( scales );
	optimizer->SetNumberOfIterations( numberOfIterations );
	optimizer->SetValueTolerance( valueTolerance );
	optimizer->SetGradientTolerance( gradientTolerance );
	optimizer->SetEpsilonFunction( epsilonFunction );

	transform->SetIdentity();
	registration->SetInitialTransformParameters( transform->GetParameters() );


	registration->SetMetric(        metric        );
	registration->SetOptimizer(     optimizer     );
	registration->SetTransform(     transform     );
	registration->SetFixedPointSet( fixed );
	registration->SetMovingPointSet(   moving   );

	try
	{
		registration->Update();
	}
	catch( itk::ExceptionObject & e )
	{
		std::cout << e << std::endl;
		exit( EXIT_FAILURE );
	}
  std::cout << "Solution = " << transform->GetParameters() 
		<<"\nMetric = "<< (registration->GetMetric()) << std::endl;

	// Write deform mesh
	// --------------------
	typedef itk::TransformMeshFilter< QEMeshType, QEMeshType,
																									TransformType > TransformMeshFilterType;
	TransformMeshFilterType::Pointer transformFilter = TransformMeshFilterType::New();
	transformFilter->SetTransform( transform );
	transformFilter->SetInput( moving );
	transformFilter->Update();

	writer->SetFileName( argv[argcOutputMesh] );
	writer->SetInput( transformFilter->GetOutput() );
	try{
		writer->Update();
	}catch( itk::ExceptionObject &err ){
		std::cerr << "Error while writing the deformed mesh with "<<
			argv[argcOutputMesh] << std::endl;
		std::cerr << err << std::endl;
		exit(EXIT_FAILURE);
	}

	msg = std::string( "Deformed mesh wrote with success as " ) + 
		std::string(argv[argcOutputMesh]);
	finalMessage( msg );

	if( saveTransform ){
		debug( "Write transforms");
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
		itk::TransformFileWriterTemplate<float>::Pointer dWriter =
			itk::TransformFileWriterTemplate<float>::New();
		itk::TransformFileWriterTemplate<float>::Pointer iWriter =
			itk::TransformFileWriterTemplate<float>::New();
#else
		itk::TransformFileWriter::Pointer dWriter = itk::TransformFileWriter::New();
		itk::TransformFileWriter::Pointer iWriter = itk::TransformFileWriter::New();
#endif

		dWriter->SetInput(transform);
		dWriter->SetFileName(argv[argcDirectTransform]);
		try{
			dWriter->Update();
		}catch( itk::ExceptionObject &err ){
			std::cerr << "Error while writing the direct transform as "
				<< argv[argcDirectTransform]
				<< std::endl;
		}

		iWriter->SetInput(transform->GetInverseTransform());
		iWriter->SetFileName(argv[argcInverseTransform]);
		try{
			iWriter->Update();
		}catch( itk::ExceptionObject &err ){
			std::cerr << "Error while writing the inverse transform as "
				<< argv[argcInverseTransform]
				<< std::endl;
		}
		msg = std::string( "Transforms wrote as (direct) ")
				+ std::string( argv[argcDirectTransform] ) + std::string(" and (inverse) ") 
				+ std::string( argv[argcInverseTransform] );
		finalMessage( msg );

	}

	std::cout << std::endl;
	exit( EXIT_SUCCESS );
}
