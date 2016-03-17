#ifndef __PHARMA__
#define __PHARMA__

/*
 * The pharmacokinetic model library. Use it to compute the 
 * concentration of contrast agent (CA) inside MRI volume of interests.
 *
 * 2014
 * Author J Plumat, UoA
 */

class PharmaModel{
	public:
		PharmaModel();
		~PharmaModel();
		
		// Perform St/S0
		double ComputeRatio( double St, double S0 );
		// estimat the CA concentration based on St,S0, relativity and T10
		// corr is PRE/POST
		double ComputeConcentration( double St, double S0, double corr, double vol );

	private:
		double m_relaxivity, m_T10;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "pharma.cxx"
#endif

#endif
