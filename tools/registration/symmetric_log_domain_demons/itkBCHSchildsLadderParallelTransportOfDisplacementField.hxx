#ifndef __itkBCHSchildsLadderParallelTransportOfDisplacementField_txx
#define __itkBCHSchildsLadderParallelTransportOfDisplacementField_txx
#include "itkBCHSchildsLadderParallelTransportOfDisplacementField.h"

#include <itkProgressAccumulator.h>

namespace itk
{

/**
 * Default constructor.
 */
template <class TInputImage, class TOutputImage>
BCHSchildsLadderParallelTransportOfDisplacementField<TInputImage, TOutputImage>
::BCHSchildsLadderParallelTransportOfDisplacementField()
{
  // Setup the number of required inputs
  this->SetNumberOfRequiredInputs( 2 );

  // By default we shouldn't be inplace
  this->InPlaceOff();

  // Set number of apprximation terms to default value
  m_NumberOfApproximationTerms = 2;

  // Declare sub filters
  m_Adder = AdderType::New();
  m_LieBracketFilterFirstOrder = LieBracketFilterType::New();
  m_LieBracketFilterSecondOrder = LieBracketFilterType::New();
  m_MultiplierByHalf = MultiplierType::New();
  m_MultiplierByTwelfth = MultiplierType::New();

  // Multipliers can always be inplace here
  m_MultiplierByHalf->InPlaceOn();
  m_MultiplierByTwelfth->InPlaceOn();

  m_MultiplierByHalf->SetConstant( 0.5 );
  m_MultiplierByTwelfth->SetConstant( 1.0 / 12.0 );
}

/**
 * Standard PrintSelf method.
 */
template <class TInputImage, class TOutputImage>
void
BCHSchildsLadderParallelTransportOfDisplacementField<TInputImage, TOutputImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Adder: " << m_Adder << std::endl;
  os << indent << "LieBracketFilterFirstOrder: " << m_LieBracketFilterFirstOrder << std::endl;
  os << indent << "LieBracketFilterSecondOrder: " << m_LieBracketFilterSecondOrder << std::endl;
  os << indent << "MultiplierByHalf: " << m_MultiplierByHalf << std::endl;
  os << indent << "MultiplierByTwelfth: " << m_MultiplierByTwelfth << std::endl;
  os << indent << "NumberOfApproximationTerms: " << m_NumberOfApproximationTerms << std::endl;
}

/**
 * GenerateData()
 */
template <class TInputImage, class TOutputImage>
void
BCHSchildsLadderParallelTransportOfDisplacementField<TInputImage, TOutputImage>
::GenerateData()
{
  InputFieldConstPointer leftField = this->GetInput(0);
  InputFieldConstPointer rightField = this->GetInput(1);

  // Create a progress accumulator for tracking the progress of minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();

  progress->SetMiniPipelineFilter(this);

  switch( m_NumberOfApproximationTerms )
    {
    case 2:
      {
      // lf + rf + 0.5*liebracket(lf,rf)
      progress->RegisterInternalFilter(m_LieBracketFilterFirstOrder, 0.5);
      progress->RegisterInternalFilter(m_MultiplierByHalf, 0.2);
      progress->RegisterInternalFilter(m_Adder, 0.3);

      m_LieBracketFilterFirstOrder->SetInput( 0, leftField );
      m_LieBracketFilterFirstOrder->SetInput( 1, rightField );

      m_MultiplierByHalf->SetInput( m_LieBracketFilterFirstOrder->GetOutput() );
      // constant set to 0.5 in constructor

      m_Adder->SetInput( 0, m_MultiplierByHalf->GetOutput() );
      m_Adder->SetInput( 1, leftField );
      m_Adder->SetInput( 2, rightField );
#if ( ITK_VERSION_MAJOR < 3 ) || ( ITK_VERSION_MAJOR == 3 && ITK_VERSION_MINOR < 13 )
      // Work-around for http://www.itk.org/Bug/view.php?id=8672
      m_Adder->InPlaceOff();
#else
      // Adder can be inplace since the 0th input is a temp field
      m_Adder->InPlaceOn();
#endif
      break;
      }
    default:
      {
      itkExceptionMacro(<< "NumberOfApproximationTerms ("
                        << m_NumberOfApproximationTerms << ") not supported");
      }
    }

  m_Adder->GraftOutput( this->GetOutput() );
  m_Adder->Update();
  this->GraftOutput( m_Adder->GetOutput() );
}

} // end namespace itk

#endif
