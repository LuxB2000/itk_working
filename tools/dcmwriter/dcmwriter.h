/*
 * =====================================================================================
 *
 *       Filename:  dcmwriter.h
 *
 *    Description:  Write a meta image
 *									TODO: write different types of image (TIFF, DCM, etc)
 *
 *        Version:  1.0
 *        Created:  15/08/14 11:26:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

#ifndef __DCMWRITER__
#define __DCMWRITER__

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesWriter.h"
#include <itkMetaDataDictionary.h>

#include "itksys/SystemTools.hxx"

//#include <sys/stat.h> // deal with file and folders
//#include <sys/types.h>

#if ITK_VERSION_MAJOR >= 4
#include "gdcmUIDGenerator.h"
#else
#include "gdcm/src/gdcmFile.h"
#include "gdcm/src/gdcmUtil.h"
#endif


template <class TInputImage, unsigned int Dimension>
class DcmWriter //: public ImageToImageFilter<TOuputImage>
{
	public:
		typedef unsigned short														OutputPixelType;
		typedef typename itk::Image<OutputPixelType, Dimension>	OutputImageType;
		typedef typename OutputImageType::Pointer		OutputImagePointerType;
		typedef TInputImage InputImageType;
		typedef typename TInputImage::Pointer InputImagePointerType;

		typedef itk::ImageFileWriter< OutputImageType > MetaWriterType;

		DcmWriter();
		DcmWriter(const char* path);
		~DcmWriter();

		// set the output filename
		void SetFileName(const char* path){
			m_path = path;
		}

		void SetInput(InputImagePointerType image);

		// could be use to write the ouputs as a DCM (mode=2)
		// default: mha (mode=1)
		void SpecifyOutputType(int mode){
			m_mode = mode;
		}
		void Update();

	private:
		typedef itk::GDCMImageIO ImageIOType;
		typedef itk::ImageSeriesWriter< InputImageType, OutputImageType >
				SeriesWriterType;
		std::string m_path;
		InputImagePointerType m_image;
		int m_mode;
		typename SeriesWriterType::DictionaryArrayType GenerateDico();



};

#ifndef ITK_MANUAL_INSTANTIATION
#include "dcmwriter.cxx"
#endif
#endif

