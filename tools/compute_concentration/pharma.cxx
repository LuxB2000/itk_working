#include "pharma.h"

PharmaModel::PharmaModel()
{
	// constant values
	// see http://imaging.bracco.com/us-en/products-and-solutions/contrast-media/multihance/relaxivity
	m_relaxivity = 5.9; // [1/(mMol*sec)]
	m_T10 = 2.359; // sec
}

PharmaModel::~PharmaModel()
{
}

double
PharmaModel::ComputeRatio( double S0, double St )
{
	return (St/S0);
}

double
PharmaModel::ComputeConcentration( double S0, double St, double cor, double vol )
{
	double c = ( (ComputeRatio(S0, St)*cor) ); // - 1 ) / ( m_relaxivity * m_T10 );
	//c = c / vol;
		return std::max(0.0,c); // if neg, return 0

}

