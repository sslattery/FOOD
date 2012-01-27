//---------------------------------------------------------------------------//
// \file TensorField_Def.hpp
// \author Stuart Slattery
// \brief Tensor field definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_DEF_HPP
#define FOOD_TENSORFIELD_DEF_HPP

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
				  const int entity_type,
				  const int entity_topology,
				  const int coord_type,
				  RCP_TensorTemplate tensor_template,
				  RCP_Unit unit,
				  const std::string &name )
    : d_comm(comm)
    , d_dofs(0)
    , d_domain(domain)
    , d_dfunckernel(dfunckernel)
    , d_entity_type(entity_type)
    , d_entity_topology(entity_topology)
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
void TensorField<Scalar>::attachToTagData( iBase_TagHandle dof_tag,
					   ErrorCode &error )
{
    d_dof_tag = dof_tag;

    error = 0;

    int num_tensor_component = 
	d_tensor_template->getNumComponents();

    int num_domain_entity = 0;
    iMesh_getNumOfTopo( d_domain->getMesh(),
			d_domain->getMeshSet(),
			d_entity_topology,
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

    iBase_EntityHandle *dof_entities = 0;
    int entities_allocated = num_domain_entity;
    int entities_size = 0;
    iMesh_getEntities( d_domain->getMesh(),
		       d_domain->getMeshSet(),
		       d_entity_type,
		       d_entity_topology,
		       &dof_entities,
		       &entities_allocated,
		       &entities_size,
		       &error );
    assert( iBase_SUCCESS == error );
    assert( num_domain_entity == entities_size );

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

    iMesh_createTag( d_domain->getMesh(),
		     &d_name[0],
		     TypeTraits<Scalar>::tag_size*num_tensor_component,
		     TypeTraits<Scalar>::tag_type,
		     &d_dof_tag,
		     &error,
		     (int) d_name.size() );
    assert( iBase_SUCCESS == error );
		     
    int num_domain_entity = 0;
    iMesh_getNumOfTopo( d_domain->getMesh(),
			d_domain->getMeshSet(),
			d_entity_topology,
			&num_domain_entity,
			&error );
    assert( iBase_SUCCESS == error );

    int dof_size = num_tensor_component*num_domain_entity;
    assert( (int) dof_array.size() == dof_size );
    d_dofs.clear();

    iBase_EntityHandle *dof_entities = 0;
    int entities_allocated = num_domain_entity;
    int entities_size = 0;
    iMesh_getEntities( d_domain->getMesh(),
		       d_domain->getMeshSet(),
		       d_entity_type,
		       d_entity_topology,
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
 * coordinates in a particular entity. 
 */
template<class Scalar>
void TensorField<Scalar>::evaluateDF( const iBase_EntityHandle entity,
				      const MDArray &coords,
				      const int is_param,
				      MDArray &dfunc_values )
{
    ErrorCode error = 0;

    // 1) Get the entity nodes and their coordinates.
    iBase_EntityHandle *element_nodes = 0;
    int adj_entity_handles_allocated = 8;
    int adj_entity_handles_size = 0;
    iMesh_getEntAdj( d_domain->getMesh(),
		     entity,
		     iBase_VERTEX,
		     &element_nodes,
		     &adj_entity_handles_allocated,
		     &adj_entity_handles_size,
		     &error );
    assert( iBase_SUCCESS == error );

    int coords_allocated = adj_entity_handles_size*3;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( d_domain->getMesh(),
			   element_nodes,
			   adj_entity_handles_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = adj_entity_handles_size;
    cell_node_dimensions[2] = coords.dimension(1);
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    // 2) Obtain pre-images of the set of evaluation points in the reference frame.
    MDArray reference_points( coords.dimension(0), coords.dimension(1) );
    Intrepid::CellTools<Scalar>::mapToReferenceFrame( reference_points,
						      coords,
						      cell_nodes,
						      *d_dfunckernel->getCellTopology(),
						      0 );

    // 3) Evaluate the basis at the pre-image set in the reference frame.
    MDArray basis_eval( d_dfunckernel->getBasis()->getCardinality(),
			coords.dimension(0) );
    d_dfunckernel->evaluateDF( basis_eval, reference_points );

    // 4) Transform evaluated basis values to the physical frame.
    MDArray transformed_eval( 1, 
			      d_dfunckernel->getBasis()->getCardinality(),
			      coords.dimension(0) );
    Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Scalar,MDArray,MDArray>( 
	transformed_eval, basis_eval );

    // 5) Evaluate the field using tensor components (the DOF for this entity).
    MDArray interpolated_vals( 1, d_dfunckernel->getBasis()->getCardinality() );
    Intrepid::FunctionSpaceTools::evaluate<Scalar,MDArray,MDArray>( 
	interpolated_vals, 
	getEntDF( entity, error), 
	transformed_eval );
    assert( iBase_SUCCESS == error );

    // 6) Integrate the evaluated basis function to give the values at the
    //    requested coordinates.
    for ( int m = 0; m < coords.dimension(0); ++m )
    {
	dfunc_values(m) = 0.0;
	for ( int n = 0; n < d_dfunckernel->getBasis()->getCardinality(); ++n )
	{
	    dfunc_values(m) += interpolated_vals(m,n);
	}
    }
    
    free( element_nodes );
    free( coord_array );
}

/*!
 * \brief Evaluate the degrees of freedom of this field at a set of
 * coordinates in a particular entity. 
 */
template<class Scalar>
void TensorField<Scalar>::evaluateGradDF( const iBase_EntityHandle entity,
					  const MDArray &coords,
					  const int is_param,
					  MDArray &dfunc_values )
{
    ErrorCode error = 0;

    // 1) Get the entity nodes and their coordinates.
    iBase_EntityHandle *element_nodes = 0;
    int adj_entity_handles_allocated = 8;
    int adj_entity_handles_size = 0;
    iMesh_getEntAdj( d_domain->getMesh(),
		     entity,
		     iBase_VERTEX,
		     &element_nodes,
		     &adj_entity_handles_allocated,
		     &adj_entity_handles_size,
		     &error );
    assert( iBase_SUCCESS == error );

    int coords_allocated = adj_entity_handles_size*3;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( d_domain->getMesh(),
			   element_nodes,
			   adj_entity_handles_size,
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = adj_entity_handles_size;
    cell_node_dimensions[2] = coords.dimension(1);
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    // 2) Obtain pre-images of the set of evaluation points in the reference frame.
    MDArray reference_points( coords.dimension(0), coords.dimension(1) );
    Intrepid::CellTools<Scalar>::mapToReferenceFrame( reference_points,
						      coords,
						      cell_nodes,
						      *d_dfunckernel->getCellTopology(),
						      0 );

    // 3) Evaluate the gradient of the basis at the pre-image set in the
    //    reference frame. 
    MDArray basis_eval( d_dfunckernel->getBasis()->getCardinality(),
			coords.dimension(0) );
    d_dfunckernel->evaluateGradDF( basis_eval, reference_points );

    // 4) Transform evaluated basis value graidents to the physical frame.
    MDArray transformed_eval( 1, 
			      d_dfunckernel->getBasis()->getCardinality(),
			      coords.dimension(0) );
    Intrepid::FunctionSpaceTools::HGRADtransformGRAD<Scalar,MDArray,MDArray>( 
	transformed_eval, basis_eval );

    // 5) Evaluate the gradient of the field using tensor components (the DOF
    // for this entity).
    MDArray interpolated_vals( 1, d_dfunckernel->getBasis()->getCardinality() );
    Intrepid::FunctionSpaceTools::evaluate<Scalar,MDArray,MDArray>( 
	interpolated_vals, 
	getEntDF( entity, error), 
	transformed_eval );
    assert( iBase_SUCCESS == error );

    // 6) Integrate the evaluated basis function gradients to give the values
    //    at the requested coordinates.
    for ( int m = 0; m < coords.dimension(0); ++m )
    {
	dfunc_values(m) = 0.0;
	for ( int n = 0; n < d_dfunckernel->getBasis()->getCardinality(); ++n )
	{
	    dfunc_values(m) += interpolated_vals(m,n);
	}
    }
    
    free( element_nodes );
    free( coord_array );
}

/*! 
 * \brief Get degrees of freedom for a particular entity in the domain.
 */
template<class Scalar>
typename TensorField<Scalar>::MDArray
TensorField<Scalar>::getEntDF( iBase_EntityHandle entity,
			       ErrorCode &error ) const
{
    error = 0;

    Teuchos::ArrayRCP<Scalar> entity_dofs( 
	d_tensor_template->getNumComponents() );

    int tag_values_allocated = entity_dofs.size()*sizeof(Scalar);
    int tag_values_size = 0;

    iMesh_getArrData( d_domain->getMesh(),
		      &entity,
		      1,
		      d_dof_tag,
		      &entity_dofs,
		      &tag_values_allocated,
		      &tag_values_size,
		      &error );
    assert( iBase_SUCCESS == error );
    assert( tag_values_allocated == tag_values_size );
    
    Teuchos::Tuple<int,2> array_dimensions;
    coeffs_dimensions[0] = 1;
    coeffs_dimensions[1] = d_tensor_template->getNumComponents;
    MDArray dof_array( Teuchos::Array<Scalar>(array_dimensions), entity_dofs );
    
    return dof_array;
}

/*! 
 * \brief Get degrees of freedom for an array of entities in the
 * domain. Returned implicitly interleaved.
 */
template<class Scalar>
typename TensorField<Scalar>::MDArray
TensorField<Scalar>::getEntArrDF( iBase_EntityHandle *entities, 
				  int num_entities,
				  ErrorCode &error ) const
{
    error = 0;

    Teuchos::ArrayRCP<Scalar> entities_dofs( 
	num_entities*d_tensor_template->getNumComponents() );

    int tag_values_allocated = entities_dofs.size()*sizeof(Scalar);
    int tag_values_size = 0;

    iMesh_getArrData( d_domain->getMesh(),
		      entities,
		      num_entities,
		      d_dof_tag,
		      &entities_dofs,
		      &tag_values_allocated,
		      &tag_values_size,
		      &error );
    assert( iBase_SUCCESS == error );
    assert( tag_values_allocated == tag_values_size );

    Teuchos::Tuple<int,2> array_dimensions;
    coeffs_dimensions[0] = num_entities;
    coeffs_dimensions[1] = d_tensor_template->getNumComponents;
    MDArray dof_array( Teuchos::Array<Scalar>(array_dimensions), entities_dofs );
    
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

