
#ifndef __SIMPLECOMMANDITERATIONUPDATE__
#define  __SIMPLECOMMANDITERATIONUPDATE__

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
    std::cout << optimizer->GetValue() << "] {";
    //std::cout << optimizer->GetCurrentPosition() << std::endl;
		vnl_matrix<double> p(3, 3);
		p[0][0] = (double) optimizer->GetCurrentPosition()[0];
    p[0][1] = (double) optimizer->GetCurrentPosition()[1];
    p[0][2] = (double) optimizer->GetCurrentPosition()[2];
		p[1][0] = (double) optimizer->GetCurrentPosition()[3];
		p[1][1] = (double) optimizer->GetCurrentPosition()[4];
		p[1][2] = (double) optimizer->GetCurrentPosition()[5];
		p[2][0] = (double) optimizer->GetCurrentPosition()[6];
		p[2][1] = (double) optimizer->GetCurrentPosition()[7];
		p[2][2] = (double) optimizer->GetCurrentPosition()[8];
		vnl_svd<double> svd(p);
		vnl_matrix<double> r(2, 2);
		r = svd.U() * vnl_transpose(svd.V());
		double angle = vcl_asin(r[1][0]);
		std::cout << "A:" << angle * 180.0 / vnl_math::pi;
		std::cout<<" T:["<<optimizer->GetCurrentPosition()[9]<<",";
		std::cout<<optimizer->GetCurrentPosition()[10]<<",";
		std::cout<<optimizer->GetCurrentPosition()[11]<<"]}"<<std::endl;
  }
};



#endif
