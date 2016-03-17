/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLocalJacobianTransform.h,v $
  Language:  C++
  Date:      $Date: 2003/01/13 04:49:23 $
  Version:   $Revision: 1.25 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkLocalJacobianTransform_h
#define __itkLocalJacobianTransform_h

#include <iostream>
#include <stdio.h>
#include <time.h>

#include "itkTransform.h"
#include "itkExceptionObject.h"

namespace itk
{
  template <class TScalarType=double,
    unsigned int NInputDimensions  = 3, 
    unsigned int NOutputDimensions = 3>
    
    class LocalJacobianTransform : 
    
    public Transform< TScalarType, 
    NInputDimensions, 
    NOutputDimensions >
    {
      public:

      /** Standard class typedefs. */
      typedef LocalJacobianTransform Self;
      typedef Transform<TScalarType,
                        NInputDimensions,
                        NOutputDimensions>             Superclass;
      typedef SmartPointer<Self>                       Pointer;
      typedef SmartPointer<const Self>                 ConstPointer;

      /** Superclass typenames*/
      typedef typename Superclass::ParametersType ParametersType;
      typedef typename Superclass::OutputPointType OutputPointType;
      typedef typename Superclass::InputPointType InputPointType;
      typedef typename Superclass::JacobianType JacobianType;
      
      /** New macro for creation of through the object factory.*/
      itkNewMacro( Self );
      
      /** Run-time type information (and related methods). */
      itkTypeMacro( LocalJacobianTransform, Transform );

      /** Gives the non-null colunms in the jacobian matrice. */
      std::vector<unsigned int> GetNonNullJacobianIndices(){
	return m_NonNullJacobianIndices;
      }
      
      /** Apply transform to a point and chek if 
	  the input point is inside the transformation domain. */
      virtual OutputPointType TransformPoint(const InputPointType  & point , bool & valid) const
      { valid = true;
	return Superclass::TransformPoint(point);
      }
      
      /** Get the jacobian and the non-null column indices. */
      virtual const JacobianType & GetJacobian(const InputPointType  &point,
					       std::vector<unsigned int> & nonNullJacobianIndices) const
      {
	return Superclass::GetJacobian(point);
      }
      

      /** Get the jacobian. */
      virtual const JacobianType & GetJacobian(const InputPointType  &point) const
      {
	return Superclass::GetJacobian(point);
      }

      protected:
      // Default constructor
      LocalJacobianTransform(){};
      // Constructor with space dimens. and number of params.
      LocalJacobianTransform
      (unsigned int dimension,unsigned int numberOfParameters)
      :Superclass(dimension,numberOfParameters){};
      ~LocalJacobianTransform(){};
      
      mutable std::vector<unsigned int> m_NonNullJacobianIndices;
      
      private:
      LocalJacobianTransform(const Self&); //purposely not implemented
      void operator=(const Self&); //purposely not implemented
      
    }; //class LocalJacobianTransform
  
}//END NAMESPACE

#endif /* __itkLocalJacobianTransform_h */
