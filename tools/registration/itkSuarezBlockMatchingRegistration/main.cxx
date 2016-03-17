
#include "itkSuarezBlockMatchingRegistration.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkNanWarpImageFilter.h"
#include "itkNumericTraits.h"

#include <cstdlib>

int main( int argc, char ** argv )
{

  ////////////////////////////////////////////////////////// PARSING

  if ( argc < 7 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[ 0 ]
              << " fixedImageFile movingImage"
              << " varianceRegularization varianceSimilarity"
              << " levelsNotToCompute"
              << " outputImageFile"
              << std::endl;
    return -1;
    }

  ////////////////////////////////////////////////////////// INPUT FILES

  typedef unsigned char PixelType;
  typedef itk::Image<PixelType, 2> ScalarImageType;

  typedef itk::Vector<float, 2> VectorType;
  typedef itk::Image<VectorType, 2> VectorImageType;
  

  typedef itk::ImageFileReader< ScalarImageType > ReaderType;
  ReaderType::Pointer fixedReader = ReaderType::New();
  fixedReader->SetFileName( argv[ 1 ] );
  try
    {
    fixedReader->Update();
    }
  catch ( itk::ExceptionObject &err )
    {
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return -1;
    }

  ReaderType::Pointer movingReader = ReaderType::New();
  movingReader->SetFileName( argv[ 2 ] );
  try
    {
    movingReader->Update();
    }
  catch ( itk::ExceptionObject &err )
    {
    std::cout << "ExceptionObject caught !" << std::endl; 
    std::cout << err << std::endl; 
    return -1;
    }

  typedef itk::ImageFileWriter<ScalarImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[ 6 ] );

  float regularization, similarity;
  regularization = atof( argv[ 3 ] );
  similarity = atof( argv[ 4 ] );

  unsigned levelsNotToCompute;
  levelsNotToCompute = atoi( argv[ 5 ] );

  ////////////////////////////////////////////////////////// PIPELINE///////////////////

  typedef itk::SuarezBlockMatchingRegistration<ScalarImageType,VectorImageType> FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetFixedImage( fixedReader->GetOutput() );
  filter->SetMovingImage( movingReader->GetOutput() );
  filter->SetVarianceRegularization( regularization );
  filter->SetVarianceSimilarity( similarity );
  filter->SetNumberOfLevelsNotToCompute( levelsNotToCompute );

  typedef itk::LinearInterpolateImageFunction<ScalarImageType, double>
    LinearInterpolatorType; // TOD: can't take 'double' off with gcc
  LinearInterpolatorType::Pointer linearInterpolator =
    LinearInterpolatorType::New();

  typedef itk::NanWarpImageFilter
    <ScalarImageType, ScalarImageType, VectorImageType> WarpType;
  WarpType::Pointer warp = WarpType::New();
  
  VectorImageType::Pointer deformation = filter->GetOutput();

  warp->SetInterpolator( linearInterpolator );
  warp->SetEdgePaddingValue(
    itk::NumericTraits<WarpType::PixelType>::quiet_NaN() );
  warp->SetOutputSpacing( deformation->GetSpacing() );
  warp->SetOutputOrigin( deformation->GetOrigin() );
  warp->SetDeformationField( deformation );
  warp->SetInput( movingReader->GetOutput() );

  writer->SetInput( warp->GetOutput() );
  writer->Update();

  return 0;

}
