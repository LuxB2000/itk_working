/*
 * =====================================================================================
 *
 *       Filename:  main.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/01/15 11:02:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 * original code: http://www.insight-journal.org/browse/publication/154, Tom Vercauteren, INRIA & Mauna Kea Technologies
 */

#include <itkCommand.h>
#include "itkDiffeomorphicDemonsRegistrationFilter.h"
#include <itkDisplacementFieldJacobianDeterminantFilter.h>
#include "itkFastSymmetricForcesDemonsRegistrationFilter.h"
#include "itkGridForwardWarpImageFilter.h"
#include <itkHistogramMatchingImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMultiResolutionPDEDeformableRegistration.h>
#include <itkTransformFileReader.h>
//#include <itkTransformToDisplacementFieldSource.h>
#include "itkVectorCentralDifferenceImageFunction.h"
#include <itkVectorLinearInterpolateNearestNeighborExtrapolateImageFunction.h>
#include "itkWarpHarmonicEnergyCalculator.h"
#include <itkWarpImageFilter.h>

#include "../../dcmwriter/dcmwriter.h"
#include "../../dcmreader/dcmreader.h"

//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
template <class TPixel=float, unsigned int VImageDimension=3>
class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate                         Self;
  typedef  itk::Command                                   Superclass;
  typedef  itk::SmartPointer<Self>                        Pointer;

  typedef itk::Image< TPixel, VImageDimension >           InternalImageType;
  typedef itk::Vector< TPixel, VImageDimension >          VectorPixelType;
  typedef itk::Image<  VectorPixelType, VImageDimension > DisplacementFieldType;

  typedef itk::DiffeomorphicDemonsRegistrationFilter<
    InternalImageType,
    InternalImageType,
    DisplacementFieldType>                                DiffeomorphicDemonsRegistrationFilterType;

  typedef itk::FastSymmetricForcesDemonsRegistrationFilter<
     InternalImageType,
     InternalImageType,
     DisplacementFieldType>                                FastSymmetricForcesDemonsRegistrationFilterType;

  typedef itk::MultiResolutionPDEDeformableRegistration<
     InternalImageType, InternalImageType,
     DisplacementFieldType, TPixel >                       MultiResRegistrationFilterType;

  typedef itk::DisplacementFieldJacobianDeterminantFilter<
     DisplacementFieldType, TPixel>                        JacobianFilterType;
  
  typedef itk::MinimumMaximumImageCalculator<
     InternalImageType>                                   MinMaxFilterType;

  typedef itk::WarpHarmonicEnergyCalculator<
     DisplacementFieldType>                                HarmonicEnergyCalculatorType;

  typedef itk::VectorCentralDifferenceImageFunction<
     DisplacementFieldType>                                WarpGradientCalculatorType;

  typedef typename WarpGradientCalculatorType::OutputType WarpGradientType;
  
  itkNewMacro( Self );

  void SetTrueField(const DisplacementFieldType * truefield)
    {
    m_TrueField = truefield;

    m_TrueWarpGradientCalculator = WarpGradientCalculatorType::New();
    m_TrueWarpGradientCalculator->SetInputImage( m_TrueField );
    
    m_CompWarpGradientCalculator =  WarpGradientCalculatorType::New();
    }
  
  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
      return;
      }

    typename DisplacementFieldType::ConstPointer deffield = 0;
    unsigned int iter = -1;
    double metricbefore = -1.0;
    
    if ( const DiffeomorphicDemonsRegistrationFilterType * dfilter = 
         dynamic_cast< const DiffeomorphicDemonsRegistrationFilterType * >( object ) )
      {
      iter = dfilter->GetElapsedIterations() - 1;
      metricbefore = dfilter->GetMetric();
      deffield = const_cast<DiffeomorphicDemonsRegistrationFilterType *>(
        dfilter)->GetDisplacementField();
      }
    else if ( const FastSymmetricForcesDemonsRegistrationFilterType * ffilter = 
              dynamic_cast< const FastSymmetricForcesDemonsRegistrationFilterType * >( object ) )
      {
      iter = ffilter->GetElapsedIterations() - 1;
      metricbefore = ffilter->GetMetric();
      deffield = const_cast<FastSymmetricForcesDemonsRegistrationFilterType *>(
        ffilter)->GetDisplacementField();
      }
    else if ( const MultiResRegistrationFilterType * multiresfilter = 
              dynamic_cast< const MultiResRegistrationFilterType * >( object ) )
      {
      std::cout<<"Finished Multi-resolution iteration :"<<multiresfilter->GetCurrentLevel()-1<<std::endl;
      std::cout<<"=============================="<<std::endl<<std::endl;
      }
    else
      {
      return;
      }

    if (deffield)
      {
      std::cout<<iter<<": MSE "<<metricbefore<<" - ";

      double fieldDist = -1.0;
      double fieldGradDist = -1.0;
      double tmp;
      if (m_TrueField)
        {
        typedef itk::ImageRegionConstIteratorWithIndex<DisplacementFieldType>
           FieldIteratorType;
        FieldIteratorType currIter(
           deffield, deffield->GetLargestPossibleRegion() );
        FieldIteratorType trueIter(
           m_TrueField, deffield->GetLargestPossibleRegion() );
        
        m_CompWarpGradientCalculator->SetInputImage( deffield );
        
        fieldDist = 0.0;
        fieldGradDist = 0.0;
        for ( currIter.GoToBegin(), trueIter.GoToBegin();
              ! currIter.IsAtEnd(); ++currIter, ++trueIter )
          {
          fieldDist += (currIter.Value() - trueIter.Value()).GetSquaredNorm();
          
          // No need to add Id matrix here as we do a substraction
          tmp = (
             ( m_CompWarpGradientCalculator->EvaluateAtIndex(currIter.GetIndex())
               -m_TrueWarpGradientCalculator->EvaluateAtIndex(trueIter.GetIndex())
                ).GetVnlMatrix() ).frobenius_norm();
          fieldGradDist += tmp*tmp;
          }
        fieldDist = sqrt( fieldDist/ (double)(
                             deffield->GetLargestPossibleRegion().GetNumberOfPixels()) );
        fieldGradDist = sqrt( fieldGradDist/ (double)(
                                 deffield->GetLargestPossibleRegion().GetNumberOfPixels()) );
        
        std::cout<<"d(.,true) "<<fieldDist<<" - ";
        std::cout<<"d(.,Jac(true)) "<<fieldGradDist<<" - ";
        }
        
      m_HarmonicEnergyCalculator->SetImage( deffield );
      m_HarmonicEnergyCalculator->Compute();
      const double harmonicEnergy
         = m_HarmonicEnergyCalculator->GetHarmonicEnergy();
      std::cout<<"harmo. "<<harmonicEnergy<<" - ";
      
      
      m_JacobianFilter->SetInput( deffield );
      m_JacobianFilter->UpdateLargestPossibleRegion();
      
      
      const unsigned int numPix = m_JacobianFilter->
         GetOutput()->GetLargestPossibleRegion().GetNumberOfPixels();
      
      TPixel* pix_start = m_JacobianFilter->GetOutput()->GetBufferPointer();
      TPixel* pix_end = pix_start + numPix;
      
      TPixel* jac_ptr;
      
      // Get percentage of det(Jac) below 0
      unsigned int jacBelowZero(0u);
      for (jac_ptr=pix_start; jac_ptr!=pix_end; ++jac_ptr)
        {
        if ( *jac_ptr<=0.0 ) ++jacBelowZero;
        }
      const double jacBelowZeroPrc = static_cast<double>(jacBelowZero)
         / static_cast<double>(numPix);
      
      
      // Get min an max jac
      const double minJac = *(std::min_element (pix_start, pix_end));
      const double maxJac = *(std::max_element (pix_start, pix_end));
      
      // Get some quantiles
      // We don't need the jacobian image
      // we can modify/sort it in place
      jac_ptr = pix_start + static_cast<unsigned int>(0.002*numPix);
      std::nth_element(pix_start, jac_ptr, pix_end);
      const double Q002 = *jac_ptr;
      
      jac_ptr = pix_start + static_cast<unsigned int>(0.01*numPix);
      std::nth_element(pix_start, jac_ptr, pix_end);
      const double Q01 = *jac_ptr;
      
      jac_ptr = pix_start + static_cast<unsigned int>(0.99*numPix);
      std::nth_element(pix_start, jac_ptr, pix_end);
      const double Q99 = *jac_ptr;
      
      jac_ptr = pix_start + static_cast<unsigned int>(0.998*numPix);
      std::nth_element(pix_start, jac_ptr, pix_end);
      const double Q998 = *jac_ptr;
      
      
      std::cout<<"max|Jac| "<<maxJac<<" - "
               <<"min|Jac| "<<minJac<<" - "
               <<"ratio(|Jac|<=0) "<<jacBelowZeroPrc<<std::endl;  

      if (this->m_Fid.is_open())
        {
        if (! m_headerwritten)
          {
          this->m_Fid<<"Iteration"
                     <<", MSE before"
                     <<", Harmonic energy"
                     <<", min|Jac|"
                     <<", 0.2% |Jac|"
                     <<", 01% |Jac|"
                     <<", 99% |Jac|"
                     <<", 99.8% |Jac|"
                     <<", max|Jac|"
                     <<", ratio(|Jac|<=0)";
          
          if (m_TrueField)
            {
            this->m_Fid<<", dist(warp,true warp)"
                       <<", dist(Jac,true Jac)";
            }
          
          this->m_Fid<<std::endl;
          
          m_headerwritten = true;
          }
        
        this->m_Fid<<iter
                   <<", "<<metricbefore
                   <<", "<<harmonicEnergy
                   <<", "<<minJac
                   <<", "<<Q002
                   <<", "<<Q01
                   <<", "<<Q99
                   <<", "<<Q998
                   <<", "<<maxJac
                   <<", "<<jacBelowZeroPrc;
        
        if (m_TrueField)
          {
          this->m_Fid<<", "<<fieldDist
                     <<", "<<fieldGradDist;
          }
        
        this->m_Fid<<std::endl;
        }
      }
    }
  
protected:   
  CommandIterationUpdate() :
     m_Fid( "metricvalues.csv" ),
     m_headerwritten(false)
    {
    m_JacobianFilter = JacobianFilterType::New();
    m_JacobianFilter->SetUseImageSpacing( true );
    m_JacobianFilter->ReleaseDataFlagOn();
    
    m_Minmaxfilter = MinMaxFilterType::New();
    
    m_HarmonicEnergyCalculator = HarmonicEnergyCalculatorType::New();
    
    m_TrueField = 0;
    m_TrueWarpGradientCalculator = 0;
    m_CompWarpGradientCalculator = 0;
    };

  ~CommandIterationUpdate()
    {
    this->m_Fid.close();
    }

private:
  std::ofstream m_Fid;
  bool m_headerwritten;
  typename JacobianFilterType::Pointer m_JacobianFilter;
  typename MinMaxFilterType::Pointer m_Minmaxfilter;
  typename HarmonicEnergyCalculatorType::Pointer m_HarmonicEnergyCalculator;
  typename DisplacementFieldType::ConstPointer m_TrueField;
  typename WarpGradientCalculatorType::Pointer m_TrueWarpGradientCalculator;
  typename WarpGradientCalculatorType::Pointer m_CompWarpGradientCalculator;
};

// ========================================
//  Main Function
// ========================================

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
	const unsigned int argcFixedImageName = 1, argcMovingImageName = 2, 
				argcOutputImageName = 3, argcOutputDeformationField = 4;
	
	const unsigned int Dimension = 3;
	
// Declare the types of the images (float or double only)
  typedef float                                  PixelType;
  typedef itk::Image< PixelType, Dimension >     ImageType;

  typedef itk::Vector< PixelType, Dimension >    VectorPixelType;
  typedef typename itk::Image
     < VectorPixelType, Dimension >              DisplacementFieldType;


  // Images we use
  typename ImageType::Pointer fixedImage = 0;
  typename ImageType::Pointer movingImage = 0;
  typename DisplacementFieldType::Pointer inputDefField = 0;


  // Set up the file readers
	typedef DcmReader< ImageType, Dimension > FixedImageReaderType;
	typedef DcmReader< ImageType, Dimension > MovingImageReaderType;
  typedef itk::ImageFileReader< DisplacementFieldType > FieldReaderType;
  typedef itk::TransformFileReader                     TransformReaderType;

	FixedImageReaderType fixedImageReader = FixedImageReaderType();
	MovingImageReaderType movingImageReader = MovingImageReaderType();
	fixedImageReader.SetFileName( argv[argcFixedImageName] );
	movingImageReader.SetFileName( argv[argcMovingImageName] );

	try{
		fixedImageReader.Update();
		fixedImage = fixedImageReader.GetOutput();
		movingImageReader.Update();
		movingImage = movingImageReader.GetOutput();
	}
	catch( itk::ExceptionObject & err ) {
			exit(EXIT_FAILURE);
	}

	// match intensities
	typedef itk::HistogramMatchingImageFilter
       <ImageType, ImageType> MatchingFilterType;
  typename MatchingFilterType::Pointer matcher = MatchingFilterType::New();
	matcher->SetInput( fixedImageReader.GetOutput() );
	matcher->SetReferenceImage( movingImageReader.GetOutput() );
	matcher->SetNumberOfHistogramLevels( 2048 );
  matcher->SetNumberOfMatchPoints( 1024 );
  matcher->ThresholdAtMeanIntensityOn();
	// Update the matcher
  try
  {
    matcher->Update();
		movingImage = matcher->GetOutput();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cout << "Could not match the input images." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
  }

	// Set up the demons filter output
  typename DisplacementFieldType::Pointer defField = 0;

	typedef typename itk::PDEDeformableRegistrationFilter
     < ImageType, ImageType, DisplacementFieldType>   BaseRegistrationFilterType;
  typename BaseRegistrationFilterType::Pointer filter;

	// params
	// TODO: use actual interesting values
	// Maximum length of an update vector (0: no restriction) - default: 2
	double maxUpdateLength = 2;
	// Type of gradient used for computing the demons force (0 is symmetrized,
	// 1 is fixed image, 2 is warped moving image, 3 is mapped moving image) - 
	// default: 0
	unsigned int gradientType = 3;
	// Smoothing sigma for the deformation field at each iteration - default:3
	double sigma = 2.5;
	// Smoothing sigma for the update field at each iteration - default: 0
	double sigmaIter = 0.05;
	// number of iteration
	unsigned int numOfLevels = 4;
	std::vector<unsigned int> numOfIterations(numOfLevels,0);
	numOfIterations[0] = 100;
	numOfIterations[1] = 80;
	numOfIterations[2] = 60;
	numOfIterations[3] = 50;

	if( true ){
		// diffeomorphic demons
		// s <- s o exp(u) (Diffeomorphic demons)
    typedef typename itk::DiffeomorphicDemonsRegistrationFilter
       < ImageType, ImageType, DisplacementFieldType>
       ActualRegistrationFilterType;
    typedef typename ActualRegistrationFilterType::GradientType GradientType;
    
    typename ActualRegistrationFilterType::Pointer actualfilter
       = ActualRegistrationFilterType::New();
    
    actualfilter->SetMaximumUpdateStepLength( maxUpdateLength );
    actualfilter->SetUseGradientType(
       static_cast<GradientType>(gradientType) );
    filter = actualfilter;
	}else{
		// s <- s + u (ITK basic implementation)
    typedef typename itk::FastSymmetricForcesDemonsRegistrationFilter
       < ImageType, ImageType, DisplacementFieldType>
       ActualRegistrationFilterType;
    typedef typename ActualRegistrationFilterType::GradientType GradientType;
    typename ActualRegistrationFilterType::Pointer actualfilter
       = ActualRegistrationFilterType::New();
    actualfilter->SetMaximumUpdateStepLength( maxUpdateLength );
    actualfilter->SetUseGradientType(
       static_cast<GradientType>(gradientType) );
    filter = actualfilter;
	}

	filter->SmoothDisplacementFieldOn();
	filter->SetStandardDeviations( sigma );
	filter->SmoothUpdateFieldOn();
  filter->SetUpdateFieldStandardDeviations( sigmaIter );

	// Create the Command observer and register it with the registration filter.
  typename CommandIterationUpdate<PixelType, Dimension>::Pointer observer =
       CommandIterationUpdate<PixelType, Dimension>::New();
	filter->AddObserver( itk::IterationEvent(), observer );

	typedef typename itk::MultiResolutionPDEDeformableRegistration<
     ImageType, ImageType, DisplacementFieldType, PixelType >   MultiResRegistrationFilterType;
  typename MultiResRegistrationFilterType::Pointer multires = MultiResRegistrationFilterType::New();

  typedef itk::VectorLinearInterpolateNearestNeighborExtrapolateImageFunction<
     DisplacementFieldType,double> FieldInterpolatorType;
  
  typename FieldInterpolatorType::Pointer VectorInterpolator =
     FieldInterpolatorType::New();

	#if ( ITK_VERSION_MAJOR > 3 ) || ( ITK_VERSION_MAJOR == 3 && ITK_VERSION_MINOR > 8 )
  multires->GetFieldExpander()->SetInterpolator(VectorInterpolator);
#endif

	multires->SetRegistrationFilter( filter );
  multires->SetNumberOfLevels( numOfLevels );
  
  multires->SetNumberOfIterations( &numOfIterations[0] );

  multires->SetFixedImage( fixedImage );
  multires->SetMovingImage( movingImage );

	typename CommandIterationUpdate<PixelType, Dimension>::Pointer multiresobserver =
       CommandIterationUpdate<PixelType, Dimension>::New();
    multires->AddObserver( itk::IterationEvent(), multiresobserver );


// Compute the deformation field
  try
    {
    multires->UpdateLargestPossibleRegion();
    }
  catch( itk::ExceptionObject& err )
    {
    std::cout << "Unexpected error." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
    }


  // The outputs
  defField = multires->GetOutput();
  defField->DisconnectPipeline();

	// warp the result
  typedef itk::WarpImageFilter
     < ImageType, ImageType, DisplacementFieldType >  WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  warper->SetInput( movingImageReader.GetOutput() );
  warper->SetOutputSpacing( fixedImage->GetSpacing() );
  warper->SetOutputOrigin( fixedImage->GetOrigin() );
  warper->SetOutputDirection( fixedImage->GetDirection() );
  warper->SetDisplacementField( defField );

  
  // Write warped image out to file
  typedef unsigned short                                OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;
  typedef itk::CastImageFilter
     < ImageType, OutputImageType >                CastFilterType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  
  typename WriterType::Pointer      writer =  WriterType::New();
  typename CastFilterType::Pointer  caster =  CastFilterType::New();
  writer->SetFileName( argv[argcOutputImageName] );
  caster->SetInput( warper->GetOutput() );
  writer->SetInput( caster->GetOutput()   );
  //writer->SetUseCompression( true );

  try
  {
    writer->Update();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cout << "Unexpected error." << std::endl;
    std::cout << err << std::endl;
    exit( EXIT_FAILURE );
  }

	if( argc >= argcOutputDeformationField + 1 ){
		// Write the deformation field as an image of vectors.
    // Note that the file format used for writing the deformation field must be
    // capable of representing multiple components per pixel. This is the case
    // for the MetaImage and VTK file formats for example.
    typedef itk::ImageFileWriter< DisplacementFieldType > FieldWriterType;
    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetFileName( argv[argcOutputDeformationField] );
    fieldWriter->SetInput( defField );
    fieldWriter->SetUseCompression( true );
		try
      {
      fieldWriter->Update();
      }
		catch( itk::ExceptionObject& err )
      {
      std::cout << "Unexpected error." << std::endl;
      std::cout << err << std::endl;
      exit( EXIT_FAILURE );
		}
	}



}
