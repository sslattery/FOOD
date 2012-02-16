//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file TensorField_Def.hpp
 * \author Stuart Slattery
 * \brief Tensor field definition.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_DEF_HPP
#define FOOD_TENSORFIELD_DEF_HPP

#include <vector>

#include "TopologyTools.hpp"

#include <Teuchos_CommHelpers.hpp>

#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
TensorField<Scalar>::TensorField( RCP_Communicator comm,
				  RCP_Domain domain,
				  RCP_DFuncKernel dfunckernel,
				  const int coord_type,
				  RCP_TensorTemplate tensor_template,
				  RCP_Unit unit,
				  const std::string &name )
    : d_comm(comm)
    , d_dofs(0)
    , d_dof_map(0)
    , d_domain(domain)
    , d_dfunckernel(dfunckernel)
    , d_coord_type(coord_type)
    , d_tensor_template(tensor_template)
    , d_unit(unit)
    , d_name(name)
    , d_dof_tag(0)
{ /* ... */ }

/*!
 * \brief Destructor.
 */
template<class Scalar>
TensorField<Scalar>::~TensorField()
{ /* ... */ }

/*! 
 * \brief Attach this field to tag data.
 *
 * Here we set a pointer to the tag data from the specified entity topology
 * and define the Tpetra map of that vector.
 */
template<class Scalar>
void TensorField<Scalar>::attachToTagData( TagHandle dof_tag,
					   ErrorCode &error )
{
    d_dof_tag = dof_tag;

    error = 0;

    int num_tensor_component = 
	d_tensor_template->getNumComponents();

    int num_domain_entity = 0;
    iMesh_getNumOfTopo( d_domain->getMesh(),
			d_domain->getMeshSet(),
			d_dfunckernel->getDFEntityTopology(),
			&num_domain_entity,
			&error );
    assert( iBase_SUCCESS == error );

    int tag_size = 0;
    iMesh_getTagSizeValues( d_domain->getMesh(),
			    dof_tag,
			    &tag_size,
			    &error );
    assert( iBase_SUCCESS == error );
    assert( tag_size == num_tensor_component*TypeTraits<Scalar>::tag_size );

    int tag_type = 0;
    iMesh_getTagType( d_domain->getMesh(),
		      dof_tag,
		      &tag_type,
		      &error );
    assert( iBase_SUCCESS == error );
    assert( tag_type == (int) TypeTraits<Scalar>::tag_type );

    EntityHandle *dof_entities = 0;
    int entities_allocated = num_domain_entity;
    int entities_size = 0;
    iMesh_getEntities( d_domain->getMesh(),
		       d_domain->getMeshSet(),
		       d_dfunckernel->getDFEntityType(),
		       d_dfunckernel->getDFEntityTopology(),
		       &dof_entities,
		       &entities_allocated,
		       &entities_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int dof_size = 
	num_tensor_component*num_domain_entity*TypeTraits<Scalar>::tag_size;
    Teuchos::ArrayRCP<Scalar> dof_array(dof_size);

    int tag_value_allocated = 
	num_tensor_component*num_domain_entity*sizeof(Scalar);
    int tag_value_size = 0;
    iMesh_getArrData( d_domain->getMesh(),
		      dof_entities,
		      entities_size,
		      dof_tag,
		      &dof_array,
		      &tag_value_allocated,
		      &tag_value_size,
		      &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,2> dof_dimensions;
    dof_dimensions[0] = num_domain_entity;
    dof_dimensions[1] = num_tensor_component;
    
    d_dofs = MDArray( Teuchos::Array<int>(dof_dimensions), dof_array );

    mapDF();

    free( dof_entities );
}

/*!
 * \brief Attach this field to array data and tag the mesh.
 */
template<class Scalar>
void TensorField<Scalar>::attachToArrayData(
    Teuchos::ArrayRCP<Scalar> dof_array, 
    int storage_order,
    ErrorCode &error )
{
    error = 0;

    int num_tensor_component = 
	d_tensor_template->getNumComponents();

    if ( d_dof_tag == 0 )
    {
	iMesh_createTag( d_domain->getMesh(),
			 &d_name[0],
			 TypeTraits<Scalar>::tag_size*num_tensor_component,
			 TypeTraits<Scalar>::tag_type,
			 &d_dof_tag,
			 &error,
			 (int) d_name.size() );
	assert( iBase_SUCCESS == error );
    }
		     
    int num_domain_entity = 0;
    iMesh_getNumOfTopo( d_domain->getMesh(),
			d_domain->getMeshSet(),
			d_dfunckernel->getDFEntityTopology(),
			&num_domain_entity,
			&error );
    assert( iBase_SUCCESS == error );

    int dof_size = num_tensor_component*num_domain_entity;
    assert( (int) dof_array.size() == dof_size );
    d_dofs.clear();

    EntityHandle *dof_entities = 0;
    int entities_allocated = num_domain_entity;
    int entities_size = 0;
    iMesh_getEntities( d_domain->getMesh(),
		       d_domain->getMeshSet(),
		       d_dfunckernel->getDFEntityType(),
		       d_dfunckernel->getDFEntityTopology(),
		       &dof_entities,
		       &entities_allocated,
		       &entities_size,
		       &error );
    assert( iBase_SUCCESS == error );
    assert( num_domain_entity == entities_size );

    if ( iBase_INTERLEAVED == storage_order )
    {
	Teuchos::Tuple<int,2> dof_dimensions;
	dof_dimensions[0] = num_domain_entity;
	dof_dimensions[1] = num_tensor_component;
	d_dofs = MDArray( Teuchos::Array<int>(dof_dimensions), dof_array );

	int tag_values_size = 
	    num_tensor_component*entities_size*sizeof(Scalar);
	iMesh_setArrData( d_domain->getMesh(),
			  dof_entities,
			  entities_size,
			  d_dof_tag,
			  &(d_dofs.getData())[0],
			  tag_values_size,
			  &error );
	assert( iBase_SUCCESS == error );
    }

    mapDF();

    free( dof_entities );
}

/*!
 * \brief Evaluate the degrees of freedom of this field at a set of
 * coordinates in a particular entity. MDArray(C,P,component).
 */
template<class Scalar>
void TensorField<Scalar>::evaluateDF( const EntityHandle entity,
				      const MDArray &coords,
				      const int is_param,
				      MDArray &dfunc_values )
{
    ErrorCode error = 0;

    // 1) Get the entity nodes and their coordinates.
    EntityHandle *element_nodes = 0;
    int element_nodes_allocated = 0;
    int element_nodes_size = 0;
    iMesh_getEntAdj( d_domain->getMesh(),
		     entity,
		     iBase_VERTEX,
		     &element_nodes,
		     &element_nodes_allocated,
		     &element_nodes_size,
		     &error );
    assert( iBase_SUCCESS == error );

    TopologyTools::MBCN2Shards( element_nodes, 
				element_nodes_size,
				d_dfunckernel->getEvalTopology() );

    int coords_allocated = element_nodes_size*3;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( d_domain->getMesh(),
			   element_nodes,
			   element_nodes_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = element_nodes_size;
    cell_node_dimensions[2] = coords.dimension(1);
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    // 2) Obtain pre-images of the set of evaluation points in the reference
    // frame. 
    MDArray reference_points( coords.dimension(0), coords.dimension(1) );
    Intrepid::CellTools<double>::mapToReferenceFrame( 
	reference_points,
	coords,
	cell_nodes,
	*d_dfunckernel->getCellTopology(),
	0 );

    // 3) Evaluate the basis at the pre-image set in the reference frame.
    MDArray basis_eval( d_dfunckernel->getBasisCardinality(),
			coords.dimension(0) );
    d_dfunckernel->evaluateValueBasis( basis_eval, reference_points );

    // 4) Transform evaluated basis values to the physical frame.
    MDArray transformed_eval( 1, 
			      d_dfunckernel->getBasisCardinality(),
			      coords.dimension(0) );
    d_dfunckernel->transformValue( transformed_eval, basis_eval );

    // 5) Evaluate the field using tensor components (the DOF for this
    // entity).
    MDArray entity_dofs = getEntDF( entity, error );
    assert( iBase_SUCCESS == error );

    MDArray component_values( 1, coords.dimension(0) );
    MDArray component_coeffs( 1, d_dfunckernel->getBasisCardinality() );
    for ( int n = 0; n < (int) d_tensor_template->getNumComponents(); ++n )
    {
	for ( int p = 0; p < coords.dimension(0); ++p )
	{
	    component_values(0,p) = 0.0;
	}

	for ( int m = 0; m < (int) d_dfunckernel->getBasisCardinality(); ++m )
	{
	    component_coeffs(0,m) = entity_dofs(0,m,n);
	}

	Intrepid::FunctionSpaceTools::evaluate<double>( component_values,
							component_coeffs, 
							transformed_eval );

	for ( int p = 0; p < coords.dimension(0); ++p )
	{
	    dfunc_values(0,p,n) = component_values(0,p);
	}
    }

    free( element_nodes );
    free( coord_array );
}

/*!
 * \brief Evaluate the gradient of the degrees of freedom of this field at a
 * set of coordinates in a particular entity. MDArray(C,P,component,spacedim).
 */
template<class Scalar>
void TensorField<Scalar>::evaluateGradDF( const EntityHandle entity,
					  const MDArray &coords,
					  const int is_param,
					  MDArray &dfunc_values )
{
    assert( FOOD_HGRAD == d_dfunckernel->getBasisFunctionSpace() );
    ErrorCode error = 0;

    // 1) Get the entity nodes and their coordinates.
    EntityHandle *element_nodes = 0;
    int element_nodes_allocated = 0;
    int element_nodes_size = 0;
    iMesh_getEntAdj( d_domain->getMesh(),
		     entity,
		     iBase_VERTEX,
		     &element_nodes,
		     &element_nodes_allocated,
		     &element_nodes_size,
		     &error );
    assert( iBase_SUCCESS == error );

    TopologyTools::MBCN2Shards( element_nodes, 
				element_nodes_size,
				d_dfunckernel->getEvalTopology() );

    int coords_allocated = element_nodes_size*3;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( d_domain->getMesh(),
			   element_nodes,
			   element_nodes_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = element_nodes_size;
    cell_node_dimensions[2] = coords.dimension(1);
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    // 2) Obtain pre-images of the set of evaluation points in the reference
    // frame. 
    MDArray reference_points( coords.dimension(0), coords.dimension(1) );
    Intrepid::CellTools<double>::mapToReferenceFrame( 
	reference_points,
	coords,
	cell_nodes,
	*d_dfunckernel->getCellTopology(),
	0 );

    // 3) Evaluate the basis at the pre-image set in the reference frame.
    MDArray basis_eval( d_dfunckernel->getBasisCardinality(),
			coords.dimension(0),
			coords.dimension(1) );
    d_dfunckernel->evaluateGradBasis( basis_eval, reference_points );

    // 4) Transform evaluated basis values to the physical frame.
    MDArray jacobian( 1, 
		      coords.dimension(0),
		      coords.dimension(1),
		      coords.dimension(1) );
    Intrepid::CellTools<double>::setJacobian( 
	jacobian, 
	reference_points,
	cell_nodes,
	*d_dfunckernel->getCellTopology() );

    MDArray jacobian_inv( 1, 
			  coords.dimension(0),
			  coords.dimension(1),
			  coords.dimension(1) );
    Intrepid::CellTools<double>::setJacobianInv( jacobian_inv, jacobian );

    MDArray transformed_eval( 1, 
			      d_dfunckernel->getBasisCardinality(),
			      coords.dimension(0),
			      coords.dimension(1) );
    Intrepid::FunctionSpaceTools::HGRADtransformGRAD<double>( 
	transformed_eval, jacobian_inv, basis_eval );

    // 5) Evaluate the field using tensor components (the DOF for this
    // entity).
    MDArray entity_dofs = getEntDF( entity, error );
    assert( iBase_SUCCESS == error );

    MDArray component_values( 1, coords.dimension(0), coords.dimension(1) );
    MDArray component_coeffs( 1, d_dfunckernel->getBasisCardinality() );
    for ( int n = 0; n < (int) d_tensor_template->getNumComponents(); ++n )
    {
	for ( int p = 0; p < coords.dimension(0); ++p )
	{
	    for ( int d = 0; d < coords.dimension(1); ++d )
	    {
		component_values(0,p,d) = 0.0;
	    }
	}

	for ( int m = 0; m < (int) d_dfunckernel->getBasisCardinality(); ++m )
	{
	    component_coeffs(0,m) = entity_dofs(0,m,n);
	}

	Intrepid::FunctionSpaceTools::evaluate<double>( 
	    component_values,
	    component_coeffs, 
	    transformed_eval );

	for ( int p = 0; p < coords.dimension(0); ++p )
	{
	    for ( int d = 0; d < coords.dimension(1); ++d )
	    {
		dfunc_values(0,p,n,d) = component_values(0,p,d);
	    }
	}
    }

    free( element_nodes );
    free( coord_array );
}

/*! 
 * \brief Get all degrees of freedom for a particular entity in the domain.
 *  MDArray(C,F,component).
 */
template<class Scalar>
typename TensorField<Scalar>::MDArray
TensorField<Scalar>::getEntDF( EntityHandle entity,
			       ErrorCode &error ) const
{
    error = 0;

    int dof_size = d_dfunckernel->getBasisCardinality()*
		   d_tensor_template->getNumComponents();

    Teuchos::ArrayRCP<Scalar> entity_dofs( dof_size );

    int tag_values_allocated = entity_dofs.size()*sizeof(Scalar);
    int tag_values_size = 0;

    EntityHandle *dof_entities = 0;
    int dof_entities_allocated = 0;
    int dof_entities_size = 0;
    iMesh_getEntAdj( d_domain->getMesh(),
		     entity,
		     d_dfunckernel->getDFEntityTopology(),
		     &dof_entities,
		     &dof_entities_allocated,
		     &dof_entities_size,
		     &error );
    assert( iBase_SUCCESS == error );

    iMesh_getArrData( d_domain->getMesh(),
		      dof_entities,
		      dof_entities_size,
		      d_dof_tag,
		      &entity_dofs,
		      &tag_values_allocated,
		      &tag_values_size,
		      &error );
    assert( iBase_SUCCESS == error );
    
    Teuchos::Tuple<int,3> array_dimensions;
    array_dimensions[0] = 1;
    array_dimensions[1] = d_dfunckernel->getBasisCardinality();
    array_dimensions[2] = d_tensor_template->getNumComponents();
    MDArray dof_array( Teuchos::Array<int>(array_dimensions), entity_dofs );
    
    return dof_array;
}

/*!
 * \brief Map the degrees of freedom with globally unique ID's.
 */
template<class Scalar>
void TensorField<Scalar>::mapDF()
{
    int myRank = d_comm->getRank();
    int local_size = d_dofs.size();
    int global_size = 0;
    Teuchos::reduceAll<int>( *d_comm,
			     Teuchos::REDUCE_MAX,
			     int(1),
			     &local_size,
			     &global_size );
    int offset = 0;
    std::vector<OrdinalType> dof_ordinals( d_dofs.size() );
    std::vector<OrdinalType>::iterator dof_ordinal_iterator;
    for ( dof_ordinal_iterator = dof_ordinals.begin();
	  dof_ordinal_iterator != dof_ordinals.end();
	  ++dof_ordinal_iterator, ++offset )
    {
	*dof_ordinal_iterator = global_size*myRank + offset;
    }
    Teuchos::ArrayView<const OrdinalType> dof_ordinals_view(dof_ordinals);
    d_dof_map = 
	Tpetra::createNonContigMap<OrdinalType>( dof_ordinals_view, d_comm );
}

} // end namespace FOOD

#endif // end FOOD_TENSORFIELD_DEF_HPP

//---------------------------------------------------------------------------//
// end TensorField_Def.hpp
//---------------------------------------------------------------------------//

