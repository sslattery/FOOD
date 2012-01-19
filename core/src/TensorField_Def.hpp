//---------------------------------------------------------------------------//
// \file TensorField_Def.hpp
// \author Stuart Slattery
// \brief Tensor field definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_DEF_HPP
#define FOOD_TENSORFIELD_DEF_HPP

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class ScalarType>
TensorField<ScalarType>::TensorField( RCP_Communicator comm,
				      RCP_Domain domain,
				      int entity_type,
				      int entity_topology,
				      int coord_type,
				      RCP_TensorTemplate tensor_template,
				      RCP_Unit unit,
				      const std::string &name )
    : d_comm(comm)
    , d_domain(domain)
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
template<class ScalarType>
TensorField<ScalarType>::~TensorField()
{ /* ... */ }

/*! 
 * \brief Attach this field to tag data.
 *
 * Here we copy the tag data from the specified entity topology into a local
 * vector and define the Tpetra map of that vector.
 */
template<class ScalarType>
void TensorField<ScalarType>::attachToTagData( iBase_TagHandle dof_tag,
					       ErrorCode &error )
{
    d_dof_tag = dof_tag;

    error = 0;

    int num_tensor_component = 
	d_tensor_template->getTensorTemplateNumComponents();

    int num_domain_entity = 0;
    iMesh_getNumOfTopo( d_domain->getDomainMesh(),
			d_domain->getDomainMeshSet(),
			d_entity_topology,
			&num_domain_entity,
			&error );
    assert( iBase_SUCCESS == error );

    int dof_size = num_tensor_component*num_domain_entity;
    d_dofs.clear();
    d_dofs.resize(dof_size);

    int tag_size = 0;
    iMesh_getTagSizeValues( d_domain->getDomainMesh(),
			    dof_tag,
			    &tag_size,
			    &error );
    assert( iBase_SUCCESS == error );
    assert( tag_size == num_tensor_component*TypeTraits<ScalarType>::tag_size );

    int tag_type = 0;
    iMesh_getTagType( d_domain->getDomainMesh(),
		      dof_tag,
		      &tag_type,
		      &error );
    assert( iBase_SUCCESS == error );
    assert( tag_type == TypeTraits<ScalarType>::tag_type );

    iBase_EntityHandle *dof_entities = 0;
    int entities_allocated = num_domain_entity;
    int entities_size = 0;
    iMesh_getEntities( d_domain->getDomainMesh(),
		       d_domain->getDomainMeshSet(),
		       d_entity_type,
		       d_entity_topology,
		       &dof_entities,
		       &entities_allocated,
		       &entities_size,
		       &error );
    assert( iBase_SUCCESS == error );
    assert( num_domain_entity == entities_size );

    int tag_value_allocated = 
	num_tensor_component*num_domain_entity*sizeof(ScalarType);
    int tag_value_size = 0;
    iMesh_getArrData( d_domain->getDomainMesh(),
		      dof_entities,
		      entities_size,
		      dof_tag,
		      &d_dofs,
		      &tag_value_allocated,
		      &tag_value_size,
		      &error );
    assert( iBase_SUCCESS == error );

    mapDF();
}

/*!
 * \brief Attach this field to array data and tag the mesh.
 */
template<class ScalarType>
void TensorField<ScalarType>::attachToArrayData( 
    Teuchos::ArrayRCP<ScalarType> dof_array, 
    int storage_order,
    ErrorCode &error )
{
    error = 0;

    int num_tensor_component = 
	d_tensor_template->getTensorTemplateNumComponents();

    int tag_size = num_tensor_component*TypeTraits<ScalarType>::tag_size;
    iMesh_createTag( d_domain->getDomainMesh(),
		     &d_name[0],
		     tag_size,
		     TypeTraits<ScalarType>::tag_type,
		     &d_dof_tag,
		     &error,
		     (int) d_name.size() );
    assert( iBase_SUCCESS == error );
		     
    int num_domain_entity = 0;
    iMesh_getNumOfTopo( d_domain->getDomainMesh(),
			d_domain->getDomainMeshSet(),
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
    iMesh_getEntities( d_domain->getDomainMesh(),
		       d_domain->getDomainMeshSet(),
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
	d_dofs = dof_array;

	int tag_values_size = 
	    num_tensor_component*entities_size*sizeof(ScalarType);
	iMesh_setArrData( d_domain->getDomainMesh(),
			  dof_entities,
			  entities_size,
			  d_dof_tag,
			  &d_dofs[0],
			  tag_values_size,
			  &error );
	assert( iBase_SUCCESS == error );
    }

    mapDF();
}

/*! 
 * \brief Get const degrees of freedom for a particular entity in the domain.
 */
template<class ScalarType>
Teuchos::ArrayRCP<const ScalarType>
TensorField<ScalarType>::getTensorFieldConstEntDF( iBase_EntityHandle entity,
						   ErrorCode &error ) const
{
    error = 0;

    Teuchos::ArrayRCP<ScalarType> 
	entity_dofs( d_tensor_template->getTensorTemplateNumComponents() );

    int tag_values_allocated = entity_dofs.size()*sizeof(ScalarType);
    int tag_values_size = 0;

    iMesh_getArrData( d_domain->getDomainMesh(),
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
template<class ScalarType>
Teuchos::ArrayRCP<const ScalarType>
TensorField<ScalarType>::getTensorFieldConstEntArrDF( 
    iBase_EntityHandle *entities, 
    int num_entities,
    ErrorCode &error ) const
{
    error = 0;

    Teuchos::ArrayRCP<ScalarType> entities_dofs( 
	num_entities*d_tensor_template->getTensorTemplateNumComponents() );

    int tag_values_allocated = entities_dofs.size()*sizeof(ScalarType);
    int tag_values_size = 0;

    iMesh_getArrData( d_domain->getDomainMesh(),
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
 * \brief Map the degrees of freedom.
 */
template<class ScalarType>
void TensorField<ScalarType>::mapDF()
{
    int myRank = d_comm->getRank();
    int mySize = d_comm->getSize();
    int offset = 0;
    std::vector<OrdinalType> dof_ordinals( d_dofs.size() );
    std::vector<OrdinalType>::iterator dof_ordinal_iterator;
    for ( dof_ordinal_iterator = dof_ordinals.begin();
	  dof_ordinal_iterator != dof_ordinals.end();
	  ++dof_ordinal_iterator, ++offset )
    {
	*dof_ordinal_iterator = myRank*mySize + offset;
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

