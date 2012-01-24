//---------------------------------------------------------------------------//
// \file TensorField_Def.hpp
// \author Stuart Slattery
// \brief Tensor field definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_DEF_HPP
#define FOOD_TENSORFIELD_DEF_HPP

#include <Teuchos_CommHelpers.hpp>

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
void TensorField<Scalar>::evaluateDF( iBase_EntityHandle entity,
				      const MDArray &coords,
				      const int is_param,
				      MDArray &dfunc_values )
{
    if( is_param )
    {
	d_dfunckernel->evaluateDF( dfunc_values, coords );
    }
    else
    {
	
    }
}

/*! 
 * \brief Get const degrees of freedom for a particular entity in the domain.
 */
template<class Scalar>
Teuchos::ArrayRCP<const Scalar>
TensorField<Scalar>::getConstEntDF( iBase_EntityHandle entity,
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
    
    return entity_dofs;
}

/*! 
 * \brief Get const degrees of freedom for an array of entities in the
 * domain. Returned implicitly interleaved.
 */
template<class Scalar>
Teuchos::ArrayRCP<const Scalar>
TensorField<Scalar>::getConstEntArrDF( iBase_EntityHandle *entities, 
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
    
    return entities_dofs;
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

