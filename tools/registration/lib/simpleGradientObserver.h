#ifndef __SIMPLEGRADIENTOBSERVER__
#define __SIMPLEGRADIENTOBSERVER__

/*
 * Observer
 */
//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"
#include "itkRegularStepGradientDescentOptimizer.h"
template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command 
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );
protected:
  RegistrationInterfaceCommand() {};
public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
  typedef   OptimizerType *                            OptimizerPointer;

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
    // If this is the first resolution level we set the maximum step length
    // (representing the first step size) and the minimum step length
		// (representing the convergence criterion) to large values.  
		// At each subsequent resolution level, we will reduce the minimum step 
		// length by a factor of 10 in order to allow the optimizer to focus
		// on progressively smaller regions. The maximum step length is set up
		// to the current step length. In this way, when the optimizer is 
		// reinitialized at the beginning of the registration process for
    // the next level, the step length will simply start with the last value
		// used for the previous level. This will guarantee the continuity of
		// the path taken by the optimizer through the parameter space.

    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( 
                       registration->GetOptimizer() );

    std::cout << "-------------------------------------" << std::endl;
    std::cout << "MultiResolution Level : "
              << registration->GetCurrentLevel() << "\n";

  }

  void Execute(const itk::Object * , const itk::EventObject & )
    { return; }
};

#endif
