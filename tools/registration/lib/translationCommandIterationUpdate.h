/*
 * =====================================================================================
 *
 *       Filename:  translationCommandIterationUpdate.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  26/11/14 16:01:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

#ifndef __ITERATIONCOMMANDITERATIONUPDATE__
#define __ITERATIONCOMMANDITERATIONUPDATE__

#include "itkCommand.h"
//#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkGradientDescentOptimizer.h"

class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdate() {};
public:
  //typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef itk::GradientDescentOptimizer OptimizerType;
  typedef   const OptimizerType *              OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    OptimizerPointer optimizer = 
      dynamic_cast< OptimizerPointer >( object );
    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }
    std::cout <<"("<< optimizer->GetCurrentIteration() << ")_[";
    std::cout << optimizer->GetValue() << "] ";
		std::cout<<" T:["<<optimizer->GetCurrentPosition()[0]<<",";
		std::cout<<optimizer->GetCurrentPosition()[1]<<",";
		std::cout<<optimizer->GetCurrentPosition()[2]<<"]"<<std::endl;
  }
};



#endif
