#ifndef __CONCENTRATIONREADER__
#define __CONCENTRATIONREADER__

/*
 * =====================================================================================
 *
 *       Filename:  dcmwriter.h
 *
 *    Description:  Class to read a concentration volume.
 *									These volumes are characterised by a vector in each
 *									physical position.
 *									The vector length is automatically detected
 *
 *        Version:  1.0
 *        Created:  14/08/14 16:46:12
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

// standard cpp
#include <string>
#include <stdio.h>


// itk dependencies
#include "itkImage.h"
#include "itkImageFileReader.h" // read META

#include "itkCastImageFilter.h"


template <class TOutputImage, unsigned int Dimension>
class ConcentrationReader //: public ImageToImageFilter<TOuputImage>
{
	public:
		typedef TOutputImage                     OutputImageType;
		typedef typename TOutputImage::Pointer   OutputPointerType;
		typedef typename TOutputImage::PixelType OutputPixelType;

		ConcentrationReader(const char* path) ;
		ConcentrationReader() ;
		~ConcentrationReader();

		// set the file name to read
		void SetFileName( const char* path){
			m_path = path;
		}

		// use it to read
		void Update();
		// use it to get the pointer to the image
		OutputPointerType GetOutput();

	private:
		typedef double InputPixelType;
		typedef typename itk::Image<InputPixelType, Dimension>	InputImageType;

		// reader types
		typedef itk::ImageFileReader< InputImageType >  ArrayReaderType;

		std::string m_path;
		OutputPointerType m_output;
		unsigned int m_type;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "concentrationreader.cxx"
#endif

#endif
