#include "meshtools.h"

/*
 * Constructor
 */
template < class TImage, class TMesh >
MeshTools<TImage,TMesh>::MeshTools()
{
}

/*
 * Constructor
 */
template < class TImage, class TMesh >
MeshTools<TImage,TMesh>::MeshTools(ImagePointerType image, 
		typename MeshType::Pointer mesh)
{
	this->SetMesh( mesh );
	this->SetVolume( image );
}

/*
 * Destructor
 */
template < class TImage, class TMesh >
MeshTools<TImage,TMesh>::~MeshTools()
{
}

/*
 * Compute the region inside the mesh
 */
template< class TImage, class TMesh >
typename MeshTools<TImage,TMesh>::ImageRegionType 
MeshTools<TImage,TMesh>::ComputeMeshRegion()
{

	const unsigned Dimension = 3;

	const typename MeshType::BoundingBoxType* bbox = m_mesh->GetBoundingBox();
	typename ImageType::IndexType start, startAbs, end, zeros;
	zeros.Fill(0);
	typename ImageType::SizeType size;
	typename ImageType::PointType pointS, pointSAbs, pointE;
	for( unsigned int i=0; i<Dimension; i++){
		pointS[i] = bbox->GetBounds()[i*2];
		pointE[i] = bbox->GetBounds()[1+i*2];
	}

	std::cout << bbox->GetBounds() << std::endl;
	std::cout << "\npointS: " << pointS << "\npointE "<< pointE << std::endl;

	for( unsigned int i=0; i<Dimension; i++ ){
		pointSAbs[i] = abs(pointS[i]);
	//	pointE[i] = abs(pointE[i]);
	}
	m_image->TransformPhysicalPointToIndex(pointS,start); 
	m_image->TransformPhysicalPointToIndex(pointSAbs,startAbs); 
	m_image->TransformPhysicalPointToIndex(pointE,end ); 

	for( unsigned int i=0; i<Dimension; i++ ){
		if( end[i] > start[i] )
			size[i] = end[i] - start[i];
		else
			size[i] = start[i] - end[i];
		//start[i] = (start[i]>0) ? start[i] : 0;
	}


	typename ImageType::IndexType orig;

	std::cout << "Start " << start << " with abs " << startAbs << std::endl;
	std::cout << "End " << end << std::endl;
	std::cout << "Size " << size << std::endl;

	m_image->TransformPhysicalPointToIndex(m_image->GetOrigin(),orig);

	startAbs[2] = start[2];
	typename ImageType::RegionType m_region(start, size);

	return m_region;

}
