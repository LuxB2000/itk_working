#ifndef __VOLUME4DREADER__
#define __VOLUME4DREADER__

/*
 */

// standard cpp
#include <string>
#include <stdio.h>


// itk dependencies
#include "itkImage.h"
#include "itkImageFileReader.h" // read META

#include "itkCastImageFilter.h"


template <class TOutputImage, unsigned int Dimension>
class Volume4dReader //: public ImageToImageFilter<TOuputImage>
{
	public:
		typedef TOutputImage                     OutputImageType;
		typedef typename TOutputImage::Pointer   OutputPointerType;
		typedef typename TOutputImage::PixelType OutputPixelType;

		Volume4dReader(const char* path) ;
		Volume4dReader() ;
		~Volume4dReader();

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
		typedef itk::VariableLengthVector<InputPixelType> ArrayType;
		typedef typename itk::VectorImage<InputPixelType, Dimension>	InputImageType;

		// reader types
		typedef itk::ImageFileReader< InputImageType >  ArrayReaderType;

		std::string m_path;
		OutputPointerType m_output;
		unsigned int m_type;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "volume4dreader.cxx"
#endif

#endif
