/*
 * =====================================================================================
 *
 *       Filename:  MRIcorrections.cxx
 *
 *    Description:  Segment the phantom based on a prior information
 *
 *        Version:  1.0
 *        Created:  21/11/14 10:25:07
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

#include "MRIcorrections.h"

/*
 * constructor
 */
template<class TImage,  const unsigned int Dimension>
MRICorrections<TImage,Dimension>::MRICorrections()
{
	// default value
	m_thresh = 1000;
}

/*
 * destructor
 */
template<class TImage,  const unsigned int Dimension>
MRICorrections<TImage,Dimension>::~MRICorrections()
{
}

/*
 * Update
 */
template<class TImage,  const unsigned int Dimension>
void
MRICorrections<TImage,Dimension>::Update()
{
	typename InputImage::SizeType szPre = m_pre->GetRequestedRegion().GetSize();
	typename InputImage::SizeType szPost = m_post->GetRequestedRegion().GetSize();

	std::cout << "test: " << szPost << std::endl;
	m_vect = OutputVectorType( szPost[2], 1 );

	typename InputImage::SizeType sz;
	for( unsigned int d=0; d<Dimension; d++ ){
		sz[d] = /*(szPost[d]<szPre[d]) ? szPost[d] : */ szPre[d] ;
	}

	std::vector<long> preN = std::vector<long>(szPre[2], 0);
	std::vector<long> postN = std::vector<long>(szPost[2], 0);
	std::vector<long> noizePreN = std::vector<long>(szPre[2], 0);
	std::vector<long> noizePostN = std::vector<long>(szPost[2], 0);
	m_preMean = std::vector<double>(sz[2], 0);
	m_postMean = std::vector<double>(sz[2], 0);
	m_noizePre_mean = std::vector<double>(sz[2], 0);
	m_noizePre_std = std::vector<double>(sz[2], 0);
	m_noizePost_mean = std::vector<double>(sz[2], 0);
	m_noizePost_std =  std::vector<double>(sz[2], 0);


	std::cout << "test: " << sz << std::endl;

	// parse the post image and segment the phantom
	// assume that the pre image is register 
	int z=0;
	long postNtmp = 0, preNtmp = 0;
	typename InputImage::IndexType ind, indPost;
	typename InputImage::PointType pt;
	for( z=0, ind[2]=0; ind[2]<sz[2]; z++, ind[2]++){
		for( ind[1]=0; ind[1]<sz[1]; ind[1]++){
			for( ind[0]=0; ind[0]<sz[0]; ind[0]++){

				// is it a phantom pixel?
				if( m_pre->GetPixel( ind ) > m_thresh ){
					preN.at(z) = preN.at(z) + 1;
					m_preMean.at(z) = m_preMean.at(z) + m_pre->GetPixel( ind );
				}else{
					noizePreN.at(z) += 1;
					m_noizePre_mean.at(z) = m_noizePre_mean.at(z) +
						m_pre->GetPixel(ind);
				}

				m_pre->TransformIndexToPhysicalPoint( ind, pt);
				m_post->TransformPhysicalPointToIndex( pt, indPost);

				if( m_post->GetPixel( indPost ) > m_thresh ){
					postN.at(z) = postN.at(z) + 1;
					m_postMean.at(z) = m_postMean.at(z) + m_post->GetPixel( indPost );
				}else{
					noizePostN.at(z) += 1;
					m_noizePost_mean.at(z) = m_noizePost_mean.at(z) +
						m_post->GetPixel( indPost );
				}
			}
		}
	}

	// compute the means
	for( z=0; z<sz[2]; z++ ){
		postNtmp += postN.at(z);
		preNtmp += preN.at(z);

		if( preN.at(z) > 0 ){
			m_preMean.at(z) = m_preMean.at(z) / ((double) preN.at(z)) ;
		}else{
			m_preMean.at(z) = 1;
		}

		if( postN.at(z) > 0 ){
			m_postMean.at(z) = m_postMean.at(z) / ((double) postN.at(z));
		}else{
			m_postMean.at(z) = 1;
		}
		m_noizePre_mean.at(z) = m_noizePre_mean.at(z) / ((double) noizePreN.at(z));
		m_noizePost_mean.at(z) = m_noizePost_mean.at(z) / ((double) noizePostN.at(z));

		// compute correction vector
		//if( (m_postMean.at(z) != 1) && (m_preMean.at(z) != 1) ){
			m_vect.at(z) = m_preMean.at(z) / m_postMean.at(z);
			if( m_vect.at(z) > 10000 || m_vect.at(z) < 0 )
				m_vect.at(z) = 1;
		//}
	}

	//
	// compute standard deviations
	//
	for( z=0, ind[2]=0; ind[2]<sz[2]; z++, ind[2]++){
		for( ind[1]=0; ind[1]<sz[1]; ind[1]++){
			for( ind[0]=0; ind[0]<sz[0]; ind[0]++){
				if( m_pre->GetPixel( ind ) > m_thresh ){
					m_noizePre_std.at(z) += pow( m_noizePre_mean.at(z)-
						m_pre->GetPixel(ind), 2);
				}
				m_pre->TransformIndexToPhysicalPoint( ind, pt);
				m_post->TransformPhysicalPointToIndex( pt, indPost);
				if( m_post->GetPixel( indPost ) > m_thresh ){
					m_noizePost_std.at(z) += pow( m_noizePost_mean.at(z)-
						m_post->GetPixel(indPost), 2);
				}
			}
		}
	}
	for( z=0, ind[2]=0; ind[2]<sz[2]; z++, ind[2]++){
		m_noizePre_std.at(z) = sqrt( m_noizePre_std.at(z) / noizePreN.at(z) );
		m_noizePost_std.at(z) = sqrt( m_noizePost_std.at(z) / noizePostN.at(z) );
	}

	typename InputImage::SizeType
		postSize = m_post->GetLargestPossibleRegion().GetSize(),
		preSize = m_pre->GetLargestPossibleRegion().GetSize();
	m_ratioPost = ((double )postNtmp) /
		((double)(postSize[0]*postSize[1]*postSize[2]));
	m_ratioPre = ( (double) preNtmp) /
		((double)(preSize[0]*preSize[1]*preSize[2]));
	
}
