#ifndef __MRICORRECTIONS__
#define __MRICORRECTIONS__

/*
 * =====================================================================================
 *
 *       Filename:  MRIcorrections.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  21/11/14 10:05:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

template <class TImage, unsigned int Dimension>
class MRICorrections{
public:
	typedef TImage InputImage;
	typedef typename InputImage::Pointer InputImagePointer;
	typedef typename InputImage::PixelType PixelType;

	typedef std::vector< double > OutputVectorType;
	MRICorrections( );
	~MRICorrections();

	// inputs
	void SetThresholdValue(PixelType val){
		m_thresh = val;
	}
	void SetPreImage(InputImagePointer preI){
		m_pre = preI;
	}
	void SetPostImage( InputImagePointer postI){
		m_post = postI;
	}
	void SetThreshold( float thresh ){
		m_thresh = thresh;
	}

	// Upadete
	void Update();

	// outputs
	OutputVectorType GetCorrectionVector(){
		return m_vect;
	}

	OutputVectorType GetPreCoeff(){
		return m_preMean;
	}

	OutputVectorType GetPostCoeff(){
		return m_postMean;
	}

	double GetPostRatio(){
		return m_ratioPost;
	}

	double GetPreRatio(){
		return m_ratioPre;
	}

	OutputVectorType GetPreNoizeMean(){
		return m_noizePre_mean;
	}
	OutputVectorType GetPreNoizeStd(){
		return m_noizePre_std;
	}
	OutputVectorType GetPostNoizeMean(){
		return m_noizePost_mean;
	}
	OutputVectorType GetPostNoizeStd(){
		return m_noizePost_std;
	}

private:
	InputImagePointer m_pre, m_post;
	OutputVectorType m_vect, m_preMean, m_postMean, 
									 m_noizePre_std, m_noizePre_mean,
									 m_noizePost_std, m_noizePost_mean;
	PixelType m_thresh;
	// ratio of the image used to compute coefficients
	double m_ratioPre, m_ratioPost;

};

#ifndef ITK_MANUAL_INSTANTIATION
#include "MRIcorrections.cxx"
#endif

#endif
