#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <stdlib.h>
#include <stdio.h>
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkShrinkImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkImage.h"
#include "itkBSplineDeformableTransform.h"
#include "itkLBFGSBOptimizer.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMetaImageIO.h"
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include "itkResampleImageFilter.h"
#include "itkCommand.h"
#include "itkSamplesMeanSquaresImageToImageMetric.h"
#include "itkSamplesNormalizedCorrelationImageToImageMetric.h"
#include "itkTimeProbesCollectorBase.h"

#ifdef ITK_USE_REVIEW
#include "itkMemoryProbesCollectorBase.h"
#define itkProbesCreate()  \
  itk::TimeProbesCollectorBase chronometer; \
  itk::MemoryProbesCollectorBase memorymeter
#define itkProbesStart( text ) memorymeter.Start( text ); chronometer.Start( text )
#define itkProbesStop( text )  chronometer.Stop( text ); memorymeter.Stop( text  )
#define itkProbesReport( stream )  chronometer.Report( stream ); memorymeter.Report( stream  )
#else
#define itkProbesCreate()  \
  itk::TimeProbesCollectorBase chronometer
#define itkProbesStart( text ) chronometer.Start( text )
#define itkProbesStop( text )  chronometer.Stop( text )
#define itkProbesReport( stream )  chronometer.Report( stream )
#endif


// Class commanditerationupdate
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////:
/*
template <class TOptimizer>
class CommandIterationUpdate : public itk::Command {
public:
  typedef CommandIterationUpdate Self;
  typedef itk::Command Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  typedef const TOptimizer * OptimizerPointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() { }
  ~CommandIterationUpdate() { }
    // the non-const execute mehtod just calls the const one
  void Execute( itk::Object * caller, const itk::EventObject & event ){
    Execute( (const itk::Object *)caller, event );}
   // const execute method
  void Execute( const itk::Object * object, const itk::EventObject & event ){
    // get the optimizer that triggered the event
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( object );
    // was there really an iteration Event?
    if ( ! itk::IterationEvent().CheckEvent( &event ) ) return;
    // print iteration and value to screen
    std::cout << optimizer->GetCurrentIteration();
    //  std::cout << " = " << (double)(optimizer->GetValue())<<std::endl;
    //std::cout << " = " << (1.0 + (double)(optimizer->GetValue()))*100.0<<std::endl;
}
};*/
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  typedef itk::LBFGSBOptimizer     OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
      Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
      OptimizerPointer optimizer = 
        dynamic_cast< OptimizerPointer >( object );
      if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
        return;
        }
      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetValue() << "   ";
      std::cout << optimizer->GetInfinityNormOfProjectedGradient() << std::endl;
    }
};

// Main
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[] )
{
  if( (argc<5) ){
    std::cerr << "Wrong number of Parameters (" << argc << ")" << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " FixedImageDir MovingImageDir OutputField Samples_percentage Number_of_grid_point Max_num_Iterations" << std::endl << std::endl;
    std::cerr << "Provided parameters :" << std::endl;
	for (int i = 1; i < argc; i++){
		std::cout << "ARG " << i << " : " << argv[i] << std::endl;}
    return 1;
    }

  std::cout << std::endl << "BEGIN  OF   THE   REGISTRATION   FROM   " << argv[2] << "   TO   " << argv[1] << std::endl;
  std::cout << "WRITING   RESULTS   IN   " << argv[3] << ".  (Number of params = " << argc << ")" << std::endl;


//Definitions and instantiations
////////////////////////////////////////////

  const    unsigned int    ImageDimension = 3;
  typedef  short           PixelType;
  typedef itk::Vector< float, ImageDimension > VectorPixelType;
  typedef itk::Image< PixelType, ImageDimension >  ImageType;
  typedef double CoordinateRepType;
  const unsigned int SplineOrder = 3;
  typedef double CoordinateRepType;
  typedef itk::BSplineDeformableTransform<  CoordinateRepType,
                            ImageDimension,  SplineOrder >     DTransformType;
  typedef itk::LBFGSBOptimizer       DOptimizerType;  
  //typedef itk::SamplesMeanSquaresImageToImageMetric< ImageType, ImageType >  DMetricType;
  typedef itk::SamplesNormalizedCorrelationImageToImageMetric< ImageType, ImageType >  DMetricType;  
  typedef itk::LinearInterpolateImageFunction< ImageType, double>  InterpolatorType;
  typedef itk::ImageRegistrationMethod< ImageType, ImageType > DRegistrationType;
  typedef DTransformType::RegionType RegionType;
  typedef DTransformType::SpacingType SpacingType;
  typedef DTransformType::OriginType OriginType;
  typedef DTransformType::ParametersType     ParametersType;
  //typedef CommandIterationUpdate< DOptimizerType > DObserver;
  typedef itk::ResampleImageFilter< ImageType, ImageType > ResampleFilterType;
  typedef itk::Vector< double, ImageDimension >  DVectorType;
  typedef itk::Image< DVectorType, ImageDimension >  DeformationFieldType;
  typedef itk::ImageRegionIterator< DeformationFieldType > FieldIterator;
  typedef itk::ImageRegionIteratorWithIndex< DeformationFieldType > DeformationIteratorType;
  typedef DOptimizerType::ScalesType       OptimizerScalesType;
  typedef itk::MetaImageIO       ImageIOType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  typedef itk::Image< PixelType, 2 >    Image2DType;
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
  typedef itk::ImageFileWriter< DeformationFieldType > FieldWriterType;


// Reading Input Images
///////////////////////////////////

  ImageIOType::Pointer imageIO = ImageIOType::New(); 
  ReaderType::Pointer fixedImageReader = ReaderType::New();
  ReaderType::Pointer movingImageReader = ReaderType::New();
  
  fixedImageReader->SetImageIO( imageIO );
  fixedImageReader->SetFileName( argv[1] );
  ImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  try{
      fixedImageReader->Update();}
    catch (itk::ExceptionObject &ex){
      std::cout << ex << std::endl;
      return EXIT_FAILURE;}   
  
  movingImageReader->SetImageIO( imageIO );
  movingImageReader->SetFileName( argv[2] );
  ImageType::Pointer movingImage = movingImageReader->GetOutput();
  try{
      movingImageReader->Update();}
    catch (itk::ExceptionObject &ex){
      std::cout << ex << std::endl;
      return EXIT_FAILURE;}

////////////////////////////////////

  ImageType::RegionType fixedRegion = fixedImage->GetLargestPossibleRegion();
  ImageType::SizeType fixedImageSize = fixedRegion.GetSize();
  ImageType::IndexType fixedImageStart = fixedRegion.GetIndex();
  ImageType::RegionType movingRegion = movingImage->GetLargestPossibleRegion();
  ImageType::SizeType movingImageSize = movingRegion.GetSize();
  ImageType::IndexType movingImageStart = movingRegion.GetIndex();
  ImageType::SpacingType fixedSpacing = fixedImage->GetSpacing();
  ImageType::SpacingType movingSpacing = movingImage->GetSpacing();

  std::cout << "Fixed image size: " << fixedImageSize << std::endl;
  std::cout << "Moving image size: " << movingImageSize << std::endl;
  std::cout << "Fixed image origin: " << fixedImage->GetOrigin() << std::endl;
  std::cout << "Moving image origin: " << movingImage->GetOrigin() << std::endl;
  std::cout << "Fixed image spacing: " << fixedSpacing << std::endl;
  std::cout << "Moving image spacing: " << movingSpacing << std::endl;
  
  
// Deformable registration
//////////////////////////

    DeformationFieldType::Pointer field = DeformationFieldType::New();
    field->SetRegions( fixedRegion );
    field->SetOrigin( fixedImage->GetOrigin() );
    field->SetSpacing( fixedImage->GetSpacing() );
    field->Allocate();
    std::cout << "Initializing accumulated field... " << std::endl;
    DeformationIteratorType initit(field, field->GetBufferedRegion());

    for (initit.GoToBegin(); !initit.IsAtEnd(); ++initit) {
	VectorPixelType v;
        for (unsigned int i=0; i<ImageDimension; i++){
                v[i] = 0.0;};
        initit.Set(v);
    }

  DMetricType::Pointer         Dmetric        = DMetricType::New();
  DOptimizerType::Pointer      Doptimizer     = DOptimizerType::New();
  DRegistrationType::Pointer   Dregistration  = DRegistrationType::New();
  InterpolatorType::Pointer    interpolator   = InterpolatorType::New(); 
  DTransformType::Pointer      Dtransform     = DTransformType::New();
  ResampleFilterType::Pointer  Dresampler     = ResampleFilterType::New();   
  
  Dmetric -> ReinitializeSeed(76926294);
if(argc>5){
  Dmetric -> SetNumberOfSamples((unsigned long) ( (double)(atoi(argv[4]))/100*(double)((  ((fixedRegion.GetSize())[0]) * ((fixedRegion.GetSize())[1]) * ((fixedRegion.GetSize())[2]) ))));
  std::cout << "Samples : " << ((unsigned long) ( (double)(atoi(argv[4]))/100*(double)((  ((fixedRegion.GetSize())[0]) * ((fixedRegion.GetSize())[1]) * ((fixedRegion.GetSize())[2]) )))) << std::endl ;}
else{
    Dmetric -> SetNumberOfSamples((unsigned long) ( 10.0/100*(double)((  ((fixedRegion.GetSize())[0]) * ((fixedRegion.GetSize())[1]) * ((fixedRegion.GetSize())[2]) ))));
  std::cout << "Samples : " << ((unsigned long) ( 10.0/100*(double)((  ((fixedRegion.GetSize())[0]) * ((fixedRegion.GetSize())[1]) * ((fixedRegion.GetSize())[2]) )))) << std::endl ;}

  Dregistration->SetMetric(Dmetric);
  Dregistration->SetOptimizer(Doptimizer);
  Dregistration->SetInterpolator(interpolator);
  Dregistration->SetTransform( Dtransform );
  Dregistration->SetFixedImage(fixedImage);
  Dregistration->SetMovingImage(movingImage);

// Region choice
//////////////////////

  Dregistration->SetFixedImageRegion(fixedRegion);
  SpacingType gridSpacing = fixedImage->GetSpacing();
  SpacingType gridSpacing_in_voxels = gridSpacing;
  OriginType gridOrigin = fixedImage->GetOrigin();
  RegionType bsplineRegion;
  RegionType::SizeType   gridSizeOnImage;
  RegionType::SizeType   gridBorderSize;
  RegionType::SizeType   totalGridSize;

  int grid_size = 3;
  int nb_grid_points = 1;
  if(argc>5){
    for(unsigned int r=0; r<ImageDimension; r++){
      gridSizeOnImage[r] = (int) ( static_cast<double>(fixedImageSize[r] - 1) / (static_cast<double>(atoi(argv[5]))/static_cast<double>(fixedSpacing[r]) )) + 1;
      nb_grid_points *= gridSizeOnImage[r];
    }
  }
  else{
    gridSizeOnImage.Fill(grid_size);
  }    
  
  gridBorderSize.Fill(SplineOrder);  
  totalGridSize = gridSizeOnImage + gridBorderSize;
  bsplineRegion.SetSize( totalGridSize );

  for(unsigned int r=0; r<ImageDimension; r++){
    gridSpacing[r] *= static_cast<double>(fixedImageSize[r] - 1)  / 
                  static_cast<double>(gridSizeOnImage[r] - 1);
    gridOrigin[r]  -=  gridSpacing[r]; 
    gridSpacing_in_voxels[r] = static_cast<double>(fixedImageSize[r] - 1)  / 
                  static_cast<double>(gridSizeOnImage[r] - 1);
  }

  Dtransform->SetGridSpacing( gridSpacing );
  Dtransform->SetGridOrigin( gridOrigin );
  Dtransform->SetGridRegion( bsplineRegion );
  std::cout << "Grid size (on image) :  " << gridSizeOnImage << std::endl;
  std::cout << "Grid spacing :  " << gridSpacing_in_voxels << "(voxels)" << std::endl;
  std::cout << "Grid spacing :  " << gridSpacing << "(mm)" << std::endl;
  std::cout << "Number of grid points : " << nb_grid_points << "pts" << std::endl;

// Setting  parameters
//////////////////////////////

  const unsigned int numberOfParameters =
  Dtransform->GetNumberOfParameters();
  ParametersType parameters( numberOfParameters );

  parameters.Fill( 0.0 );

  Dtransform->SetParameters( parameters );

  Dregistration->SetInitialTransformParameters( Dtransform->GetParameters() );
  DOptimizerType::BoundSelectionType boundSelect( Dtransform->GetNumberOfParameters() );
  DOptimizerType::BoundValueType upperBound( Dtransform->GetNumberOfParameters() );
  DOptimizerType::BoundValueType lowerBound( Dtransform->GetNumberOfParameters() );

  boundSelect.Fill( 0 );
  upperBound.Fill( 0.0 );
  lowerBound.Fill( 0.0 );

  Doptimizer->SetBoundSelection( boundSelect );
  Doptimizer->SetUpperBound( upperBound );
  Doptimizer->SetLowerBound( lowerBound );

  int numit = 20;
  if(argc>6){
    numit = atoi(argv[6]);
  }

  Doptimizer->SetCostFunctionConvergenceFactor( 1e+10 );//1e+3
  Doptimizer->SetProjectedGradientTolerance( 1e-10 );//1e-7
  Doptimizer->SetMaximumNumberOfIterations( numit );
  Doptimizer->SetMaximumNumberOfEvaluations( numit*10 );
  Doptimizer->SetMaximumNumberOfCorrections( numit );

  // Create the Command observer and register it with the optimizer.
  //
  /*DObserver::Pointer Dobserver = DObserver::New();
  Doptimizer->AddObserver( itk::IterationEvent(), Dobserver );
  itk::TimeProbesCollectorBase collector;*/
  CommandIterationUpdate::Pointer Dobserver = CommandIterationUpdate::New();
  Doptimizer->AddObserver( itk::IterationEvent(), Dobserver );
  itkProbesCreate();

  std::cout << "Running Registration..." << std::endl;

  try{ 
    //collector.Start( "Registration" );
    itkProbesStart( "Registration" );
    Dregistration->StartRegistration(); 
    //collector.Stop( "Registration" );
    itkProbesStop( "Registration" );} 
  catch( itk::ExceptionObject & err ){ 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return -1;}       
  
// Write the output field
//////////////////////////////
    
    DOptimizerType::ParametersType DfinalParameters = 
                    Dregistration->GetLastTransformParameters();

    //collector.Report();
    itkProbesReport( std::cout );
    Dtransform->SetParameters( DfinalParameters );
    
    std::cout << "Computing the dense deformation field ... " << std::endl;
    DeformationFieldType::Pointer newfield = DeformationFieldType::New();
    newfield->SetRegions( fixedRegion );
    newfield->SetOrigin( fixedImage->GetOrigin() );
    newfield->SetSpacing( fixedImage->GetSpacing() );
    newfield->Allocate();
    FieldIterator fi( newfield, fixedRegion );
    fi.GoToBegin();
    DeformationFieldType::IndexType index;
    while( ! fi.IsAtEnd() ){
     index = fi.GetIndex();
     DeformationFieldType::PointType fixedPoint;
     newfield->TransformIndexToPhysicalPoint( index, fixedPoint );
     DeformationFieldType::PointType T1 = Dtransform->TransformPoint( fixedPoint );
     VectorPixelType newv;
    for (unsigned int i = 0; i < ImageDimension; i++ ){
      newv[i] = T1[i] - fixedPoint[i];}
      fi.Set( newv );
      ++fi;}

    std::cerr << "Writing deformation field in "<< argv[3] << std::endl;    
    FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    std::string fieldname = argv[3];
    fieldWriter->SetFileName( fieldname.c_str() );
    fieldWriter->SetInput( newfield );
    fieldWriter->Update();
    

  std::cout << std::endl << "END   OF   THE   REGISTRATION   FROM   " << argv[2] << "   TO   " << argv[1] << std::endl << std::endl << std::endl;

  return 0;
}
