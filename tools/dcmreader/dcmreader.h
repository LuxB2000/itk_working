#ifndef DCMREADER
#define DCMREADER

/*
 * Simple class that choose the correct way to read DCM or META
 * based on the path and the extension of the file.
 * The following extensions are took into consideration: mha, mhd for META
 * images and tif for 3D TIFF image.
 * If the path designates a folder than all the files inside the the folder
 * are assumed to be a stack of 2D slices.
 *
 * Warning: the class outputs an image as specified in the template but
 * the read images are assumed to be FLOAT if META and UNSINCHED CHAR if DCM.
 * Warning: the tiff images are assumed to be "empty of parameters" to read
 * these parameters, we are looking for a txt file with the same name (i.e. 
 * the first number in the file name) but cituated
 * m_path/../Parameters/ID-parameters.txt
 */

// standard cpp
#include <string>
#include <stdio.h>
#include <sys/stat.h> // deal with file and folders
#include <sys/types.h>


// itk dependencies
#include "itkImage.h"
#include "itkImageFileReader.h" // read META
#include "itkGDCMSeriesFileNames.h" //rad series
#include "itkImageSeriesReader.h"
#include "itkGDCMImageIO.h"
#include "itkTIFFImageIO.h"

#include "itkCastImageFilter.h"

// reading the file parameters
#include "readTiffImageParameters.h"

template <class TOutputImage, unsigned int Dimension>
class DcmReader //: public ImageToImageFilter<TOuputImage>
{
	public:
		typedef TOutputImage                     OutputImageType;
		typedef typename TOutputImage::Pointer   OutputPointerType;
		typedef typename TOutputImage::PixelType OutputPixelType;

		DcmReader(const char* path) ;
		DcmReader() ;
		~DcmReader();

		// set the file name to read
		void SetFileName( const char* path){
			m_path = path;
		}

		// use it to read
		void Update();
		// use it to get the pointer to the image
		OutputPointerType GetOutput();

		// if TIFF (and only if TIFF image) we can have access to the 
		// path to reach the parameters file. It can be used with 
		// ReadTiffImageParameters parametersReader =
		// ReadTiffImageParameters( path );
		// If not tiff image (or if the parameter file is not found)
		// it return an empty string
		std::string GetParamPath(){
			return m_param_path;
		}

	private:
		typedef unsigned short				DcmPixelType;
		typedef unsigned short				MetaPixelType;
		typedef unsigned int TIFFPixelType;
		typedef typename itk::Image<DcmPixelType,  3>		DCMImageType;
		typedef typename itk::Image<MetaPixelType, 3>		MetaImageType;
		typedef typename itk::Image<TIFFPixelType, 3>		TIFFImageType;

		// reader types
		typedef std::vector< std::string > SeriesIdContainer;
		typedef std::vector< std::string > FileNamesContainer;
		typedef itk::GDCMSeriesFileNames NamesGeneratorType;
		typedef itk::GDCMImageIO											 ImageIOType;
		typedef itk::ImageSeriesReader< DCMImageType > DCMReaderType;
		typedef itk::ImageFileReader< MetaImageType >  MetaReaderType;
		typedef itk::ImageFileReader< TIFFImageType >  TIFFReaderType;

		void m_parseName();

		std::string m_path, m_param_path;
		OutputPointerType m_output;
		unsigned int m_type;
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "dcmreader.cxx"
#endif

#endif
