//---------------------------------------------------------------------------//
// \file Domain.hpp
// \author Stuart Slattery
// \brief Domain definition
//---------------------------------------------------------------------------//

#ifndef FOOD_DOMAIN_HPP
#define FOOD_DOMAIN_HPP

#include "Types.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>

namespace FOOD
{

class Domain
{
  public:

    //@{
    //! Typedefs.
    typedef Teuchos::RCP<iMesh_Instance>               RCP_Mesh;
    typedef Teuchos::RCP<iBase_EntitySetHandle>        RCP_EntitySet;
    typedef iBase_EntityIterator*                      EntityIterator;
    typedef iBase_ErrorType                            ErrorCode;
    //@}

  private:

    // iMesh instance to which the domain set belongs.
    RCP_Mesh d_mesh;

    // Mesh set over which the field is defined.
    RCP_EntitySet d_mesh_set;

  public:

    // Constructor.
    Domain(iMesh_Instance mesh_instance, iBase_EntitySetHandle set_handle);

    // Destructor.
    ~Domain();

    //! Get the mesh instance.
    RCP_Mesh mesh() const
    { return d_mesh; }

    //! Get the mesh set.
    RCP_EntitySet meshSet() const
    { return d_mesh_set; }

    // Initialize an iterator over the specified entity type and topology.
    EntityIterator initEntIter(const int entity_type, 
			       const int entity_topology);

    // Initalize an array iterator over the specified entity type and
    // topology. 
    EntityIterator initEntArrIter(const int entity_type,
				  const int entity_topology,
				  const int array_size);
};

}

#endif // end FOOD_DOMAIN_HPP

//---------------------------------------------------------------------------//
// end Domain.hpp
//---------------------------------------------------------------------------//
