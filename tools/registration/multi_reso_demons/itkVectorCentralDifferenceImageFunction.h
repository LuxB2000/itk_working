/*=========================================================================
  *
  * Copyright Insight Software Consortium
  *
  * Licensed under the Apache License, Version 2.0 (the "License");
  * you may not use this file except in compliance with the License.
  * You may obtain a copy of the License at
  *
  * http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 * Unless required by applicable law or agreed to in writing, software
  * distributed under the License is distributed on an "AS IS" BASIS,
  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  * See the License for the specific language governing permissions and
  * limitations under the License.
  *
  *=========================================================================*/
 #ifndef __itkVectorCentralDifferenceImageFunction_h
 #define __itkVectorCentralDifferenceImageFunction_h
 
#include "itkImageFunction.h"
#include "itkMatrix.h"

namespace itk
{

template< typename TInputImage, typename TCoordRep = float >
class VectorCentralDifferenceImageFunction:
public ImageFunction< TInputImage, Matrix< double, TInputImage::PixelType::Dimension,TInputImage::ImageDimension >, TCoordRep >
	{
		public:
		typedef typename TInputImage::PixelType InputPixelType;

		itkStaticConstMacro(Dimension, unsigned int,
		InputPixelType::Dimension);

		itkStaticConstMacro(ImageDimension, unsigned int,
		TInputImage::ImageDimension);

		typedef VectorCentralDifferenceImageFunction Self;
		typedef ImageFunction< TInputImage,
		Matrix< double, itkGetStaticConstMacro(Dimension),
		itkGetStaticConstMacro(ImageDimension) >,
		TCoordRep > Superclass;
		typedef SmartPointer< Self > Pointer;
		typedef SmartPointer< const Self > ConstPointer;

		itkTypeMacro(VectorCentralDifferenceImageFunction, ImageFunction);

		itkNewMacro(Self);

		typedef TInputImage InputImageType;
		typedef typename Superclass::OutputType OutputType;
		typedef typename Superclass::IndexType IndexType;

		typedef typename Superclass::ContinuousIndexType ContinuousIndexType;
		typedef typename Superclass::PointType PointType;

		virtual OutputType EvaluateAtIndex(const IndexType & index) const;

		virtual OutputType Evaluate(const PointType & point) const
		{
			IndexType index;
			this->ConvertPointToNearestIndex(point, index);
			return this->EvaluateAtIndex(index);
		}

		virtual OutputType EvaluateAtContinuousIndex(
		const ContinuousIndexType & cindex) const
		{
			IndexType index;
			this->ConvertContinuousIndexToNearestIndex(cindex, index);
			return this->EvaluateAtIndex(index);
		}

		itkSetMacro(UseImageDirection, bool);
		itkGetConstMacro(UseImageDirection, bool);
		itkBooleanMacro(UseImageDirection);

protected:
		VectorCentralDifferenceImageFunction();
		~VectorCentralDifferenceImageFunction(){}
		void PrintSelf(std::ostream & os, Indent indent) const;

private:
		VectorCentralDifferenceImageFunction(const Self &); //purposely not
		// implemented
		void operator=(const Self &); //purposely not
		// implemented

		// flag to take or not the image direction into account
		// when computing the derivatives.
		bool m_UseImageDirection;
};
} // end namespace itk
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVectorCentralDifferenceImageFunction.hxx"
#endif

#endif

