/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  test the registration with demons - defautl example: DeformatbleRegistraion3 in ITK4
 *
 *        Version:  1.0
 *        Created:  13/10/14 16:31:29
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSymmetricForcesDemonsRegistrationFilter.h"
#include "itkHistogramMatchingImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkWarpImageFilter.h"

class CommandIterationUpdate : public itk::Command
{
	public:
		typedef CommandIterationUpdate Self;
		typedef itk::Command Superclass;
		typedef itk::SmartPointer<CommandIterationUpdate> Pointer;
		itkNewMacro( CommandIterationUpdate );
	protected:
		CommandIterationUpdate() {};
		typedef itk::Image< float, 2 > InternalImageType;
		typedef itk::Vector< float, 2 > VectorPixelType;
		typedef itk::Image< VectorPixelType, 2 > DisplacementFieldType;
		typedef itk::SymmetricForcesDemonsRegistrationFilter<
		InternalImageType,
		InternalImageType,
		DisplacementFieldType> RegistrationFilterType;
	public:
		void Execute(itk::Object *caller, const itk::EventObject & event)
		{
			Execute( (const itk::Object *)caller, event);
		}
		void Execute(const itk::Object * object, const itk::EventObject & event)
		{
			const RegistrationFilterType * filter =
				static_cast< const RegistrationFilterType * >( object );
		if( !(itk::IterationEvent().CheckEvent( &event )) ) {
			return;
		}
		std::cout << filter->GetMetric() << std::endl;
		}
};


// main
int main( int argc, char *argv[] )
{
	std::cout << "--- " << argv[0] << "---" << std::endl;
	if( argc < 4 )
	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " fixedImageFile movingImageFile ";
		std::cerr << " outputImageFile [DisplacementFieldPath] " << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "Fixed Image: " << argv[1] <<"\n"
		<< "Moving Image: " << argv[2] << "\n" 
		<< "Output Image: " << argv[3] << std::endl;

	const unsigned int Dimension = 3;
	typedef unsigned short PixelType;
	typedef itk::Image< PixelType, Dimension > FixedImageType;
	typedef itk::Image< PixelType, Dimension > MovingImageType;
	typedef itk::ImageFileReader< FixedImageType > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
	FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
	fixedImageReader->SetFileName( argv[1] );
	movingImageReader->SetFileName( argv[2] );

	typedef float InternalPixelType;
	typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
	typedef itk::CastImageFilter< FixedImageType,
	InternalImageType > FixedImageCasterType;
	typedef itk::CastImageFilter< MovingImageType,
	InternalImageType > MovingImageCasterType;
	FixedImageCasterType::Pointer fixedImageCaster = FixedImageCasterType::New();
	MovingImageCasterType::Pointer movingImageCaster
						= MovingImageCasterType::New();
	fixedImageCaster->SetInput( fixedImageReader->GetOutput() );
	movingImageCaster->SetInput( movingImageReader->GetOutput() );

	typedef itk::HistogramMatchingImageFilter<
	InternalImageType,
	InternalImageType > MatchingFilterType;
	MatchingFilterType::Pointer matcher = MatchingFilterType::New();
	matcher->SetInput( movingImageCaster->GetOutput() );
	matcher->SetReferenceImage( fixedImageCaster->GetOutput() );
	matcher->SetNumberOfHistogramLevels( 1024 );
	matcher->SetNumberOfMatchPoints( 2000 );
	matcher->ThresholdAtMeanIntensityOn();

	typedef itk::Vector< float, Dimension > VectorPixelType;
	typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldType;
	typedef itk::SymmetricForcesDemonsRegistrationFilter<
	InternalImageType,
	InternalImageType,
	DisplacementFieldType> RegistrationFilterType;
	RegistrationFilterType::Pointer filter = RegistrationFilterType::New();
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	filter->AddObserver( itk::IterationEvent(), observer );
	filter->SetFixedImage( fixedImageCaster->GetOutput() );
	filter->SetMovingImage( matcher->GetOutput() );
	filter->SetNumberOfIterations( 150 );
	filter->SetStandardDeviations( 2.0 );
	filter->Update();

	typedef itk::WarpImageFilter<
	MovingImageType,
	MovingImageType,
	DisplacementFieldType > WarperType;
	typedef itk::LinearInterpolateImageFunction<
	MovingImageType,
	double > InterpolatorType;
	WarperType::Pointer warper = WarperType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();
	warper->SetInput( movingImageReader->GetOutput() );
	warper->SetInterpolator( interpolator );
	warper->SetOutputSpacing( fixedImage->GetSpacing() );
	warper->SetOutputOrigin( fixedImage->GetOrigin() );
	warper->SetOutputDirection( fixedImage->GetDirection() );
	warper->SetDisplacementField( filter->GetOutput() );
	typedef unsigned short OutputPixelType;
	typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
	typedef itk::CastImageFilter<
	MovingImageType,
	OutputImageType > CastFilterType;
	typedef itk::ImageFileWriter< OutputImageType > WriterType;
	WriterType::Pointer writer = WriterType::New();
	CastFilterType::Pointer caster = CastFilterType::New();
	writer->SetFileName( argv[3] );
	caster->SetInput( warper->GetOutput() );
	writer->SetInput( caster->GetOutput() );
	writer->Update();
	std::cout << "Deformed image wrote as " << argv[3] << std::endl;
	if( argc > 4 ) // if a fourth line argument has been provided...
	{
		typedef itk::ImageFileWriter< DisplacementFieldType > FieldWriterType;
		FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
		fieldWriter->SetFileName( argv[4] );
		fieldWriter->SetInput( filter->GetOutput() );
		fieldWriter->Update();
		std::cout << "Deformation field wrote as " << argv[4] << std::endl;
	}
	std::cout << std::endl;
	return EXIT_SUCCESS;
}
