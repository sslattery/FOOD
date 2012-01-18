//---------------------------------------------------------------------------//
// \file TensorField_Def.hpp
// \author Stuart Slattery
// \brief Tensor field definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_TENSORFIELD_DEF_HPP
#define FOOD_TENSORFIELD_DEF_HPP

#include <valarray>

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
typename TensorField<ScalarType>::ErrorCode 
TensorField<ScalarType>::attachToTagData( iBase_TagHandle dof_tag )
{
    ErrorCode error = 0;

    int num_domain_entity = 0;
    iMesh_getNumOfTopo( d_domain->getDomainMesh(),
			d_domain->getDomainMeshSet(),
			d_entity_topology,
			&num_domain_entity,
			&error );
    assert( iBase_SUCCESS == error );

    int num_tensor_component = 
	d_tensor_template->getTensorTemplateNumComponents();

    int dof_size = num_tensor_component*num_domain_entity;
    d_dofs.clear();
    d_dofs.resize(dof_size);

    int tag_size_type = 0;
    iMesh_getTagSizeValues( d_domain->getDomainMesh(),
			    dof_tag,
			    &tag_size_type,
			    &error );
    assert( iBase_SUCCESS == error );
    assert( tag_size_type == num_tensor_component );

    int tag_size_bytes = 0;
    iMesh_getTagSizeBytes( d_domain->getDomainMesh(),
			   dof_tag,
			   &tag_size_bytes,
			   &error );
    assert( iBase_SUCCESS == error );
    assert( tag_size_bytes == (int) (num_tensor_component*sizeof(ScalarType)) );

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

    int myRank = d_comm->getRank();
    int mySize = d_comm->getSize();
    int offset = 0;
    std::vector<OrdinalType> dof_ordinals(dof_size);
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

    return error;
}

} // end namespace FOOD

#endif // end FOOD_TENSORFIELD_DEF_HPP

//---------------------------------------------------------------------------//
// end TensorField_Def.hpp
//---------------------------------------------------------------------------//

