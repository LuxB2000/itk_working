/*
 * =====================================================================================
 *
 *       Filename:  dcmwriter.h
 *
 *    Description:  Write a 4D meta image
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

#ifndef __CONCENTRATIONWRITER__
#define __CONCENTRATIONWRITER__

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"


template <class TInputImage, unsigned int Dimension>
class ConcentrationWriter //: public ImageToImageFilter<TOuputImage>
{
	public:
		typedef double														OutputPixelType;
		typedef typename itk::Image<OutputPixelType, Dimension>	OutputImageType;
		typedef typename OutputImageType::Pointer		OutputImagePointerType;
		typedef TInputImage InputImageType;
		typedef typename TInputImage::Pointer InputImagePointerType;

		typedef itk::ImageFileWriter< OutputImageType > MetaWriterType;

		ConcentrationWriter();
		ConcentrationWriter(const char* path);
		~ConcentrationWriter();

		// set the output filename
		void SetFileName(const char* path){
			m_path = path;
		}

		void SetInput(InputImagePointerType image);

		void Update();

	private:
		std::string m_path;
		InputImagePointerType m_image;
		int m_mode;



};

#ifndef ITK_MANUAL_INSTANTIATION
#include "concentrationwriter.cxx"
#endif

#endif

