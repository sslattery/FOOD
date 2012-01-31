//---------------------------------------------------------------------------//
// \file Octree.cpp
// \author Stuart Slattery
// \breif Octree declaration.
//---------------------------------------------------------------------------//

#include <cassert>
#include <algorithm>

#include "Octree.hpp"

namespace FOOD
{

/*! 
 * \brief Constructor.
 */
Octree::Octree( RCP_Domain domain, const int entity_topology )
    : d_domain(domain)
    , d_root_node( Teuchos::rcp(new OctreeNode) )
{
    d_root_node.node_set = d_domain->getMeshSet();
    d_root_node.cutting_axis = X_Axis;
    getEntSetBox( d_root_node.node_set, d_root_node.bounding_box );   
}

/*!
 * \brief Destructor.
 */
Octree::~Octree()
{ /* ... */ }

/*! 
 * \brief Build the tree for a specific entity topology.
 */
void Octree::buildTree( RCP_Node node )
{
}

/*!
 * \brief Get the bounding box of a set of entities.
 */
void Octree::getEntSetBox( const iBase_EntitySetHandle entity_set,
			   Box &bounding_box )
{
    int error = 0;

    int num_domain_vertices = 0;
    iMesh_getNumOfTopo( d_domain->getDomain()->getDomainMesh(),
			entity_set,
			iMesh_POINT,
			&num_set_vertices,
			&error );
    assert( iBase_SUCCESS == error );

    iBase_EntityHandle *set_vertices = 0;
    int set_vertices_allocated = num_set_vertices;
    int set_vertices_size = 0;
    iMesh_getEntities( d_domain->getDomain()->getDomainMesh(),
		       entity_set,
		       iBase_VERTEX,
		       iMesh_POINT,
		       &set_vertices,
		       &set_vertices_allocated,
		       &set_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int set_coords_allocated = 3*set_vertices_size;
    int set_coords_size = 0;
    Teuchos::ArrayRCP<double> set_coords( set_coords_allocated, 0.0 );
    iMesh_getVtxArrCoords( d_domain->getDomain()->getDomainMesh(),
			   entity_set,
			   iBase_BLOCKED,
			   &set_coords,
			   &set_coords_allocated,
			   &set_coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::ArrayRCP<double>::const_iterator set_x_it_begin 
	= set_coords.begin();
    Teuchos::ArrayRCP<double>::const_iterator set_x_it_end 
	= set_x_it_begin + set_coords_allocated;

    Teuchos::ArrayRCP<double>::const_iterator set_y_it_begin 
	= set_x_it_end;
    Teuchos::ArrayRCP<double>::const_iterator set_y_it_end 
	= set_y_it_begin + set_coords_allocated;
    
    Teuchos::ArrayRCP<double>::const_iterator set_z_it_begin 
	= set_y_it_end;
    Teuchos::ArrayRCP<double>::const_iterator set_z_it_end 
	= set_z_it_begin + set_coords_allocated;
    
    bounding_box[0] = *(std::min( set_x_it_begin, set_x_it_end ));
    bounding_box[1] = *(std::max( set_x_it_begin, set_x_it_end ));
    bounding_box[2] = *(std::min( set_y_it_begin, set_y_it_end ));
    bounding_box[3] = *(std::max( set_y_it_begin, set_y_it_end ));
    bounding_box[4] = *(std::min( set_z_it_begin, set_z_it_end ));
    bounding_box[5] = *(std::max( set_z_it_begin, set_z_it_end ));
}

} // end namespace food

//---------------------------------------------------------------------------//
// end Octree.cpp
//---------------------------------------------------------------------------//
