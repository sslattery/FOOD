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

#include <Teuchos_CommHelpers.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
TensorField<Scalar>::TensorField( RCP_Communicator comm,
				  RCP_Domain domain,
				  RCP_DFuncKernel dfunckernel,
				  RCP_TensorTemplate tensor_template,
				  RCP_Unit unit,
				  const std::string &name )
    : d_comm(comm)
    , d_dofs(0)
    , d_dof_map(0)
    , d_domain(domain)
    , d_dfunckernel(dfunckernel)
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
					   int &error )
{
    d_dof_tag = dof_tag;

    error = 0;

    int num_tensor_component = 
	d_tensor_template->getNumComponents();

    int num_domain_entity = 0;
    iMesh_getNumOfTopo( d_domain->getMesh(),
			d_domain->getMeshSet(),
			iMesh_POINT,
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
		       iBase_VERTEX,
		       iMesh_POINT,
		       &dof_entities,
		       &entities_allocated,
		       &entities_size,
		       &error );
    assert( iBase_SUCCESS == error );

    int dof_size = 
	num_tensor_component*num_domain_entity*TypeTraits<Scalar>::tag_size;
    d_dofs.clear();
    d_dofs.resize( dof_size );

    int tag_value_allocated = 
	num_tensor_component*num_domain_entity*sizeof(Scalar);
    int tag_value_size = 0;
    iMesh_getArrData( d_domain->getMesh(),
		      dof_entities,
		      entities_size,
		      dof_tag,
		      &d_dofs,
		      &tag_value_allocated,
		      &tag_value_size,
		      &error );
    assert( iBase_SUCCESS == error );

    mapDF();

    free( dof_entities );
}

/*!
 * \brief Attach this field to array data and tag the mesh.
 */
template<class Scalar>
void TensorField<Scalar>::attachToArrayData(
    const Teuchos::ArrayRCP<Scalar> &dof_array, 
    const int storage_order,
    int &error )
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
			iMesh_POINT,
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
		       iBase_VERTEX,
		       iMesh_POINT,
		       &dof_entities,
		       &entities_allocated,
		       &entities_size,
		       &error );
    assert( iBase_SUCCESS == error );
    assert( num_domain_entity == entities_size );

    if ( iBase_INTERLEAVED == storage_order )
    {
	d_dofs = dof_array;

	int tag_values_size = 
	    num_tensor_component*entities_size*sizeof(Scalar);
	iMesh_setArrData( d_domain->getMesh(),
			  dof_entities,
			  entities_size,
			  d_dof_tag,
			  d_dofs.get(),
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
				      const double coords[3],
				      const int is_param,
				      Teuchos::ArrayRCP<Scalar> &dfunc_values )
{
    // 1) Obtain pre-images of the set of evaluation points in the reference
    // frame. 
    double param_coords[3] = { coords[0], coords[1], coords[2] };
    if ( !is_param )
    {
	d_dfunckernel->transformPoint( param_coords,
				       coords,
				       d_domain->getMesh(),
				       entity );
    }

    // 2) Get the basis values at the pre-image set in the reference frame.
    Teuchos::ArrayRCP<Scalar> values;
    d_dfunckernel->dfuncValue( values, param_coords );

    // 3) Transform evaluated basis values to the physical frame.
    Teuchos::ArrayRCP<Scalar> transformed_values;
    d_dfunckernel->transformValue( transformed_values,
				   values,
				   param_coords,
				   d_domain->getMesh(),
				   entity );

    // 4) Evaluate the field using tensor components (the DOF for this
    // entity).
    int num_component = d_tensor_template->getNumComponents();
    int dimension = values.size()/d_dfunckernel->getCardinality();
    int dfunc_values_size = num_component*dimension;
			    
    dfunc_values.resize( dfunc_values_size );
    Teuchos::ArrayRCP<Scalar> entity_coeffs = getEntDF( entity );
    Teuchos::ArrayRCP<Scalar> component_coeffs(d_dfunckernel->getCardinality());
    Teuchos::ArrayRCP<Scalar> component_values;
    for ( int n = 0; n < num_component; ++n )
    {
	component_values = dfunc_values.persistingView( n*dimension,
							dimension );

	for ( int d = 0; d < dimension; ++d )
	{
	    for ( int m = 0; m < d_dfunckernel->getCardinality(); ++m )
	    {
		component_coeffs[m] = entity_coeffs[d + m*dimension];
	    }
	    d_dfunckernel->evaluate( component_values,
				     component_coeffs,
				     transformed_values );
	}
    }
}

/*!
 * \brief Evaluate the gradient of the degrees of freedom of this field at a
 * set of coordinates in a particular entity. 
 */
template<class Scalar>
void TensorField<Scalar>::evaluateGradDF( const iBase_EntityHandle entity,
					  const double coords[3],
					  const int is_param,
					  Teuchos::ArrayRCP<Scalar> &dfunc_values )
{
    // 1) Obtain pre-images of the set of evaluation points in the reference
    // frame. 
    double param_coords[3] = { coords[0], coords[1], coords[2] };
    if ( !is_param )
    {
	d_dfunckernel->transformPoint( param_coords,
				       coords,
				       d_domain->getMesh(),
				       entity );
    }

    // 2) Get the basis operators at the pre-image set in the reference frame.
    Teuchos::ArrayRCP<Scalar> operators;
    d_dfunckernel->dfuncOperator( operators, param_coords );

    // 3) Transform evaluated basis operators to the physical frame.
    Teuchos::ArrayRCP<Scalar> transformed_operators;
    d_dfunckernel->transformOperator( transformed_operators,
				      operators,
				      param_coords,
				      d_domain->getMesh(),
				      entity );

    // 4) Evaluate the field using tensor components (the DOF for this
    // entity).
    int dimension = operators.size()/d_dfunckernel->getCardinality();
    int dfunc_values_size = d_tensor_template->getNumComponents()*dimension;
			    
    dfunc_values.resize( dfunc_values_size );
    Teuchos::ArrayRCP<Scalar> entity_coeffs = getEntDF( entity );
    Teuchos::ArrayRCP<Scalar> component_operators;

    for ( int n = 0; n < (int) d_tensor_template->getNumComponents(); ++n )
    {
	component_operators = dfunc_values.persistingView( n*dimension, 
							   dimension );
	d_dfunckernel->evaluate( component_operators,
				 entity_coeffs,
				 transformed_operators );
    }
}

/*! 
 * \brief Get all degrees of freedom for a particular entity in the domain.
 */
template<class Scalar>
Teuchos::ArrayRCP<Scalar>
TensorField<Scalar>::getEntDF( iBase_EntityHandle entity ) const
{
    int error = 0;

    int dof_size = d_dfunckernel->getCardinality()*
		   d_tensor_template->getNumComponents();

    Teuchos::ArrayRCP<Scalar> entity_dofs( dof_size );

    iBase_EntityHandle *dof_entities = 0;
    int dof_entities_allocated = 0;
    int dof_entities_size = 0;
    iMesh_getEntAdj( d_domain->getMesh(),
		     entity,
		     iMesh_POINT,
		     &dof_entities,
		     &dof_entities_allocated,
		     &dof_entities_size,
		     &error );
    assert( iBase_SUCCESS == error );

    int tag_values_allocated = entity_dofs.size()*sizeof(Scalar);
    int tag_values_size = 0;
    iMesh_getArrData( d_domain->getMesh(),
		      dof_entities,
		      dof_entities_size,
		      d_dof_tag,
		      &entity_dofs,
		      &tag_values_allocated,
		      &tag_values_size,
		      &error );
    assert( iBase_SUCCESS == error );
    
    free( dof_entities );

    return entity_dofs;
}

/*! 
 * \brief Get all degrees of freedom for an array of entities in the domain.
 */
template<class Scalar>
Teuchos::ArrayRCP<Scalar>
TensorField<Scalar>::getEntArrDF( iBase_EntityHandle *entities,
				  int num_entities ) const
{
    int error = 0;

    int num_dof_ents = num_entities*d_dfunckernel->getCardinality();
    std::vector<iBase_EntityHandle> total_dof_entities;

    int dof_size = num_dof_ents*d_tensor_template->getNumComponents();
    Teuchos::ArrayRCP<Scalar> entity_dofs( dof_size );

    for ( int n = 0; n < num_entities; ++n )
    {
	iBase_EntityHandle *dof_entities = 0;
	int dof_entities_allocated = 0;
	int dof_entities_size = 0;
	iMesh_getEntAdj( d_domain->getMesh(),
			 entities[n],
			 iMesh_POINT,
			 &dof_entities,
			 &dof_entities_allocated,
			 &dof_entities_size,
			 &error );
	assert( iBase_SUCCESS == error );

	for ( int i = 0; i < dof_entities_size; ++i )
	{
	    total_dof_entities.push_back( dof_entities[i] );
	}

	free( dof_entities );
    }

    int tag_values_allocated = entity_dofs.size()*sizeof(Scalar);
    int tag_values_size = 0;
    iMesh_getArrData( d_domain->getMesh(),
		      total_dof_entities,
		      (int) total_dof_entities.size(),
		      d_dof_tag,
		      &entity_dofs,
		      &tag_values_allocated,
		      &tag_values_size,
		      &error );
    assert( iBase_SUCCESS == error );
    
    return entity_dofs;
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

