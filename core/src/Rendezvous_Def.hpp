//---------------------------------------------------------------------------//
// \file Rendezvous_Def.hpp
// \author Stuart Slattery
// \brief Rendezvous class definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_RENDEZVOUS_DEF_HPP
#define FOOD_RENDEZVOUS_DEF_HPP

#include <cassert>
#include <algorithm>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ArrayRCP.hpp>

#include <Tpetra_Vector.hpp>

namespace FOOD
{
 
/*!
 * \brief Constructor.
 */
template<class ScalarType>
Rendezvous<ScalarType>::Rendezvous(RCP_TensorField domain, 
				   RCP_TensorField range )
    : d_domain_primary(domain)
    , d_range_primary(range)
    , d_domain_secondary(0)
    , d_range_secondary(0)
    , d_domain_export(0)
    , d_range_export(0)
{
    assert( d_domain_primary->getTensorTemplate() == 
	    d_range_primary->getTensorTemplate() );
    assert( d_domain_primary->getComm() == 
	    d_range_primary->getComm() );
}

/*!
 * \brief Destructor.
 */
template<class ScalarType>
Rendezvous<ScalarType>::~Rendezvous()
{ /* ... */ }

/*!
 * \brief Do parallel rendezvous to generate secondary decompositions.
 */
template<class ScalarType>
void Rendezvous<ScalarType>::createSecondaryDecompositions()
{
    int error = 0;

    // 1) Compute a global axis aligned bounding box that bounds the
    // intersection of the domain and range mesh sets. 
    Teuchos::Tuple<double,6> intersection_box 
	= computeIntersectionBoundingBox();

    // 2) Create a rendezvous decomposition by performing RCB on the
    // intersecting domain and range elements.

    // 3) Create secondary decomposition fields.

    // 4) Setup export objects.
    
}

/*!
 * \brief Copy the domain degrees of freedom from the primary decomposition to
    the secondary.
*/
template<class ScalarType>
void Rendezvous<ScalarType>::domainCopyPrimaryToSecondary()
{
    Tpetra::Vector<ScalarType> primary_vector( 
	d_domain_primary->getDFMap(),
	d_domain_primary->getConstDFView() );

    Tpetra::Vector<ScalarType> secondary_vector( 
	d_domain_secondary->getDFMap() );

    secondary_vector.doExport( primary_vector, 
			       *d_domain_export, 
			       Tpetra::INSERT );
}

/*!
 * \brief Copy the range degrees of freedom from the secondary decomposition
    to the primary.
*/
template<class ScalarType>
void Rendezvous<ScalarType>::rangeCopySecondaryToPrimary()
{
    Tpetra::Vector<ScalarType> secondary_vector( 
	d_range_secondary->getDFMap(),
	d_range_secondary->getConstDFView() );

    Tpetra::Vector<ScalarType> primary_vector( 
	d_range_primary->getDFMap() );

    primary_vector.doExport( secondary_vector, 
			     *d_range_export, 
			     Tpetra::INSERT );
}

/*!
 * \brief Compute a global axis aligned bounding box that bounds the
 * intersection of the domain and range mesh sets. 
 */
template<class ScalarType>
Teuchos::Tuple<double,6> 
Rendezvous<ScalarType>::computeIntersectionBoundingBox()
{
    int error = 0;

    int num_domain_vertices = 0;
    iMesh_getNumOfTopo( d_domain_primary->getDomain()->getDomainMesh(),
			d_domain_primary->getDomain()->getDomainMeshSet(),
			iMesh_POINT,
			&num_domain_vertices,
			&error );
    assert( iBase_SUCCESS == error );

    iBase_EntityHandle *domain_vertices = 0;
    int domain_vertices_allocated = num_domain_vertices;
    int domain_vertices_size = 0;
    iMesh_getEntities( d_domain_primary->getDomain()->getDomainMesh(),
		       d_domain_primary->getDomain()->getDomainMeshSet(),
		       iBase_VERTEX,
		       iMesh_POINT,
		       &domain_vertices,
		       &domain_vertices_allocated,
		       &domain_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int domain_coords_allocated = 3*domain_vertices_size;
    int domain_coords_size = 0;
    Teuchos::ArrayRCP<double> domain_coords( domain_coords_allocated, 0.0 );
    iMesh_getVtxArrCoords( d_domain_primary->getDomain()->getDomainMesh(),
			   d_domain_primary->getDomain()->getDomainMeshSet(),
			   iBase_BLOCKED,
			   &domain_coords,
			   &domain_coords_allocated,
			   &domain_coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::ArrayRCP<double>::const_iterator domain_x_it_begin 
	= domain_coords.begin();
    Teuchos::ArrayRCP<double>::const_iterator domain_x_it_end 
	= domain_x_it_begin + domain_coords_allocated;

    Teuchos::ArrayRCP<double>::const_iterator domain_y_it_begin 
	= domain_x_it_end;
    Teuchos::ArrayRCP<double>::const_iterator domain_y_it_end 
	= domain_y_it_begin + domain_coords_allocated;
    
    Teuchos::ArrayRCP<double>::const_iterator domain_z_it_begin 
	= domain_y_it_end;
    Teuchos::ArrayRCP<double>::const_iterator domain_z_it_end 
	= domain_z_it_begin + domain_coords_allocated;
    
    Teuchos::Tuple<double,6> local_domain_box;

    local_domain_box[0] = *(std::min( domain_x_it_begin, domain_x_it_end ));
    local_domain_box[1] = *(std::max( domain_x_it_begin, domain_x_it_end ));
    local_domain_box[2] = *(std::min( domain_y_it_begin, domain_y_it_end ));
    local_domain_box[3] = *(std::max( domain_y_it_begin, domain_y_it_end ));
    local_domain_box[4] = *(std::min( domain_z_it_begin, domain_z_it_end ));
    local_domain_box[5] = *(std::max( domain_z_it_begin, domain_z_it_end ));

    Teuchos::Tuple<double,6> global_domain_box;

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MIN,
			     int(1),
			     &local_domain_box[0],
			     &global_domain_box[0] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MAX,
			     int(1),
			     &local_domain_box[1],
			     &global_domain_box[1] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MIN,
			     int(1),
			     &local_domain_box[2],
			     &global_domain_box[2] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MAX,
			     int(1),
			     &local_domain_box[3],
			     &global_domain_box[3] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MIN,
			     int(1),
			     &local_domain_box[4],
			     &global_domain_box[4] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MAX,
			     int(1),
			     &local_domain_box[5],
			     &global_domain_box[5] );

    int num_primary_range_vertices = 0;
    iMesh_getNumOfTopo( d_range_primary->getDomain()->getDomainMesh(),
			d_range_primary->getDomain()->getDomainMeshSet(),
			iMesh_POINT,
			&num_primary_range_vertices,
			&error );
    assert( iBase_SUCCESS == error );

    iBase_EntityHandle *range_vertices = 0;
    int range_vertices_allocated = num_primary_range_vertices;
    int range_vertices_size = 0;
    iMesh_getEntities( d_range_primary->getDomain()->getDomainMesh(),
		       d_range_primary->getDomain()->getDomainMeshSet(),
		       iBase_VERTEX,
		       iMesh_POINT,
		       &range_vertices,
		       &range_vertices_allocated,
		       &range_vertices_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int range_coords_allocated = 3*range_vertices_size;
    int range_coords_size = 0;
    Teuchos::ArrayRCP<double> range_coords( range_coords_allocated, 0.0 );
    iMesh_getVtxArrCoords( d_range_primary->getDomain()->getDomainMesh(),
			   d_range_primary->getDomain()->getDomainMeshSet(),
			   iBase_BLOCKED,
			   &range_coords,
			   &range_coords_allocated,
			   &range_coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::ArrayRCP<double>::const_iterator range_x_it_begin 
	= range_coords.begin();
    Teuchos::ArrayRCP<double>::const_iterator range_x_it_end 
	= range_x_it_begin + range_coords_allocated;

    Teuchos::ArrayRCP<double>::const_iterator range_y_it_begin 
	= range_x_it_end;
    Teuchos::ArrayRCP<double>::const_iterator range_y_it_end 
	= range_y_it_begin + range_coords_allocated;
    
    Teuchos::ArrayRCP<double>::const_iterator range_z_it_begin 
	= range_y_it_end;
    Teuchos::ArrayRCP<double>::const_iterator range_z_it_end 
	= range_z_it_begin + range_coords_allocated;

    Teuchos::Tuple<double,6> local_range_box;

    local_range_box[0] = *(std::min( range_x_it_begin, range_x_it_end ));
    local_range_box[1] = *(std::max( range_x_it_begin, range_x_it_end ));
    local_range_box[2] = *(std::min( range_y_it_begin, range_y_it_end ));
    local_range_box[3] = *(std::max( range_y_it_begin, range_y_it_end ));
    local_range_box[4] = *(std::min( range_z_it_begin, range_z_it_end ));
    local_range_box[5] = *(std::max( range_z_it_begin, range_z_it_end ));

    Teuchos::Tuple<double,6> global_range_box;

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MIN,
			     int(1),
			     &local_range_box[0],
			     &global_range_box[0] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MAX,
			     int(1),
			     &local_range_box[1],
			     &global_range_box[1] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MIN,
			     int(1),
			     &local_range_box[2],
			     &global_range_box[2] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MAX,
			     int(1),
			     &local_range_box[3],
			     &global_range_box[3] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MIN,
			     int(1),
			     &local_range_box[4],
			     &global_range_box[4] );

    Teuchos::reduceAll<int>( *d_domain_primary->getComm(),
			     Teuchos::REDUCE_MAX,
			     int(1),
			     &local_range_box[5],
			     &global_range_box[5] );

    return boxIntersectionTest( global_domain_box, global_range_box );
}

/*!
 * \brief Intersection test for two axis aligned boxes. Return the resulting
 * intersecting axis aligned box. 
 */
template<class ScalarType>
Teuchos::Tuple<double,6> 
Rendezvous<ScalarType>::boxIntersectionTest( Teuchos::Tuple<double,6> box_A,
					     Teuchos::Tuple<double,6> box_B )
{
    Teuchos::Tuple<double,6> intersection_box;
    
    return intersection_box;
}

} // end namespace FOOD

#endif // end FOOD_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end Rendezvous_Def.hpp
//---------------------------------------------------------------------------//
