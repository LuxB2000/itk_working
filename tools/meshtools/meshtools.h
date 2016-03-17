#ifndef __MESHTOOLS__
#define __MESHTOOLS__
/*
 * =====================================================================================
 *
 *       Filename:  meshTools.h
 *
 *    Description:  Some tools to deal with mesh
 *
 *        Version:  1.0
 *        Created:  12/08/14 09:51:39
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jerome Plumat (JP), j.plumat@auckland.ac.nz
 *        Company:  UoA, Auckand, NZ
 *
 * =====================================================================================
 */

template< class TImage, class TMesh>
class MeshTools{
	public:
		typedef TImage	ImageType;
		typedef typename TImage::Pointer	ImagePointerType;
		typedef typename ImageType::RegionType ImageRegionType;

		typedef TMesh		MeshType;
		typedef typename MeshType::Pointer MeshPointerType;

		MeshTools();
		~MeshTools();
		MeshTools( ImagePointerType image, MeshPointerType mesh);

		void SetMesh( MeshPointerType mesh){
			m_mesh=mesh;
		}

		void SetVolume( ImagePointerType image ){
			m_image = image;
		}

		ImageRegionType ComputeMeshRegion();

private:
		MeshPointerType m_mesh;
		ImagePointerType m_image;
};


#ifndef ITK_MANUAL_INSTANTIATION
#include "meshtools.cxx"
#endif


#endif
