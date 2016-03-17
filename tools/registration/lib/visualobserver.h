#ifndef __VISUALOBSERVER__
#define __VISUALOBSERVER__
/*
 * ============================================================================
 *
 *       Filename:  visualobserver.h
 *
 *    Description:  An observer that display the current result of the
 *									registration
 *
 *        Version:  1.0
 *        Created:  16/10/14 10:28:26
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * ============================================================================
 */
#include "itkExtractImageFilter.h"
#include "itkCommand.h"
#include "QuickView.h"
template <typename TFixedImage, typename TMovingImage, typename TTransform,
				 typename TRegistration, typename TOpitmizer>
class RegistrationVisualizeInterfaceCommand : public itk::Command 
{
public:
  typedef  RegistrationVisualizeInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );
	typedef TFixedImage											FixedImageType;
	typedef TMovingImage										MovingImageType;
  typedef TRegistration										RegistrationType;
  typedef RegistrationType *							RegistrationPointer;
	typedef TOpitmizer											OptimizerType;
	typedef TOpitmizer*											OptimizerPointerType;
	typedef TTransform											TransformType;

	void SetFixedImage( typename FixedImageType::Pointer f){
		m_fixed = f;
	}
	void SetMovingImage( typename MovingImageType::Pointer m){
		m_moving = m;
	}


  void Execute(itk::Object * object, const itk::EventObject & event)
	{
    // First we verify if that the event invoked is of the right type.
    // If not, we return without any further action.
    if( !(itk::IterationEvent().CheckEvent( &event )) )
    {
			return;
    }

    RegistrationPointer registration =
                            dynamic_cast<RegistrationPointer>( object );

    OptimizerPointerType optimizer = dynamic_cast< OptimizerPointerType >( 
                       registration->GetOptimizer() );

		// get the current registration parameters
		typedef typename OptimizerType::ParametersType					ParametersType;

		//
		// first, set the correct parameters values
		//
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : "
              << registration->GetCurrentLevel() << "\n";

		double maxStepLength, minStepLength;
    if ( registration->GetCurrentLevel() == 0 )
      {
				maxStepLength = 0.25; // TODO modify, initialy it was 0.5
				minStepLength = 0.01;
				/*
      optimizer->SetMaximumStepLength( 1.0 );  
      optimizer->SetMinimumStepLength( 0.01 );
			*/
      }
    else
      {
				maxStepLength = optimizer->GetMaximumStepLength() / 2;
				minStepLength = optimizer->GetMinimumStepLength() / 5.0;
				/*
      optimizer->SetMaximumStepLength( 
					optimizer->GetMaximumStepLength() / 1.75 );

      optimizer->SetMinimumStepLength( 
					optimizer->GetMinimumStepLength() / 5.0 );
					*/
      }
		std::cout << "Max/Min Step Length: " << maxStepLength 
			<< "/" << minStepLength << "\n";
		
		std::cout << std::endl;
      optimizer->SetMaximumStepLength( maxStepLength );
      optimizer->SetMinimumStepLength( minStepLength );

		//
		// display current solution
		//
		ParametersType currentSolution = registration->GetLastTransformParameters();
		typename TransformType::Pointer currentTransform = TransformType::New();
		std::cout << " ==================== " << std::endl;
		std::cout << "TEST: current solution has size: " << currentSolution.GetSize() << std::endl;
		if( currentSolution.GetSize() > 1 )
		{
			currentTransform->SetParameters( currentSolution );
			std::cout << " ==================== " << std::endl;

			typedef itk::ResampleImageFilter< 
															MovingImageType, 
															FixedImageType >			ResampleFilterType;
			typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
			resampler->SetTransform( currentTransform );
			resampler->SetInput( m_moving );

			resampler->SetSize( m_fixed->GetLargestPossibleRegion().GetSize() );
			resampler->SetOutputOrigin(  m_fixed->GetOrigin() );
			resampler->SetOutputSpacing( m_fixed->GetSpacing() );
			resampler->SetOutputDirection( m_fixed->GetDirection() );
			resampler->SetDefaultPixelValue( 0 );

			resampler->Update();

			// extract one slice
			typedef itk::Image<typename MovingImageType::PixelType, 2> SliceType;
			typedef itk::ExtractImageFilter<MovingImageType,
																				SliceType> ExtractSliceFilterType;
			typename MovingImageType::SizeType size =
				resampler->GetOutput()->GetLargestPossibleRegion().GetSize();
			size[2] = 0;
			typename MovingImageType::IndexType index = {{0,0,2}};

			typename MovingImageType::RegionType extractionRegion(index,size);
			typename ExtractSliceFilterType::Pointer extracSlice =
										ExtractSliceFilterType::New();
			extracSlice->SetInput( resampler->GetOutput() );
			extracSlice->SetExtractionRegion( extractionRegion );
			extracSlice->Update();

			typename ExtractSliceFilterType::Pointer extracFSlice =
										ExtractSliceFilterType::New();
			extracFSlice->SetInput( m_fixed );
			extracFSlice->SetExtractionRegion( extractionRegion );
			extracFSlice->Update();

			// display the slice
			//extracSlice->GetOutput()->Print( std::cout );
			std::cout << "SHOULD DISPLAY IMAGE" << std::endl;
			QuickView viewer;
			viewer.AddImage( extracSlice->GetOutput(), true);
			viewer.AddImage( extracFSlice->GetOutput(), true);
			viewer.Visualize();
		} //end if size > 1
	}//end Execute

  // Another version of the \code{Execute()} method accepting a \code{const}
  // input object is also required since this method is defined as
	// pure virtual in the base class.
	// This version simply returns without taking any action.
  void Execute(const itk::Object * , const itk::EventObject & )
    { return; }

protected:
		RegistrationVisualizeInterfaceCommand(){};


private:
		typename FixedImageType::Pointer m_fixed;
		typename MovingImageType::Pointer m_moving;

};

#endif
