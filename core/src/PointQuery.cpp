//---------------------------------------------------------------------------//
// \file PointQuery.cpp
// \author Stuart Slattery
// \brief Point query method definitions for reference elements.
//---------------------------------------------------------------------------//

#include <vector>

#include "PointQuery.hpp"
#include "CellTopologyFactory.hpp"
#include "TopologyTools.hpp"

#include <Teuchos_ENull.hpp>

#include <Intrepid_CellTools.hpp>

namespace FOOD
{

bool PointQuery::pointInRefElement( const iMesh_Instance mesh,
				    const iBase_EntityHandle entity,
				    const MDArray &coords)
{
    int error = 0;
    int topology = 0;
    iMesh_getEntTopo( mesh, entity, &topology, &error );
    assert( iBase_SUCCESS == error );

    int num_linear_nodes = TopologyTools::numLinearNodes( topology );

    iBase_EntityHandle *element_nodes = 0;
    int element_nodes_allocated = 0;
    int element_nodes_size = 0;
    iMesh_getEntAdj( mesh,
		     entity,
		     iBase_VERTEX,
		     &element_nodes,
		     &element_nodes_allocated,
		     &element_nodes_size,
		     &error );
    assert( iBase_SUCCESS == error );

    std::vector<iBase_EntityHandle> linear_nodes(num_linear_nodes);
    for ( int i = 0; i < num_linear_nodes; ++i )
    {
	linear_nodes[i] = element_nodes[i];
    }

    CellTopologyFactory topo_factory;
    RCP_CellTopology cell_topo = topo_factory.create( topology, 
						      num_linear_nodes );

    int coords_allocated = 0;
    int coords_size = 0;
    double *coord_array = 0;
    iMesh_getVtxArrCoords( mesh,
			   &linear_nodes[0],
			   (int) linear_nodes.size(),
			   iBase_INTERLEAVED,
			   &coord_array,
			   &coords_allocated,
			   &coords_size,
			   &error );
    assert( iBase_SUCCESS == error );

    Teuchos::Tuple<int,3> cell_node_dimensions;
    cell_node_dimensions[0] = 1;
    cell_node_dimensions[1] = (int) linear_nodes.size();
    cell_node_dimensions[2] = 3;
    MDArray cell_nodes( Teuchos::Array<int>(cell_node_dimensions), 
			coord_array );

    MDArray reference_points(1,3);
    Intrepid::CellTools<double>::mapToReferenceFrame( 
	reference_points,
	coords,
	cell_nodes,
	*cell_topo,
	0 );

    bool return_val = Intrepid::CellTools<double>::checkPointInclusion( 
	&reference_points[0],
	3,
	*cell_topo);

    free( element_nodes );
    free( coord_array );
    cell_topo = Teuchos::null;

    return return_val;
}

} // end namespace FOOD

//---------------------------------------------------------------------------//
// end PointQuery.cpp
//---------------------------------------------------------------------------//
