//---------------------------------------------------------------------------//
// \file PointQuery.hpp
// \author Stuart Slattery
// \brief Point query methods declaration for reference elements
//---------------------------------------------------------------------------//

#ifndef FOOD_POINTQUERY_HPP
#define FOOD_POINTQUERY_HPP

#include <cassert>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>

#include <Shards_CellTopology.hpp>

#include <Intrepid_FieldContainer.hpp>

namespace FOOD
{

namespace PointQuery
{

typedef Intrepid::FieldContainer<double>          MDArray;
typedef Teuchos::RCP<shards::CellTopology>        RCP_CellTopology;

bool pointInRefElement( const iMesh_Instance mesh,
			const iBase_EntityHandle entity,
			const MDArray &coords );
		      
} // end namespace PointQuery

} // end namespace FOOD

#endif // end FOOD_POINTQUERY_HPP

//---------------------------------------------------------------------------//
// end PointQuery.hpp
//---------------------------------------------------------------------------//
