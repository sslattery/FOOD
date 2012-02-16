//---------------------------------------------------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
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

//@{
//! Typedefs.
typedef Intrepid::FieldContainer<double>          MDArray;
typedef Teuchos::RCP<shards::CellTopology>        RCP_CellTopology;
typedef iBase_EntityHandle                        EntityHandle;
//@}

// Point in volume query.
bool pointInRefElement( const iMesh_Instance mesh,
			const EntityHandle entity,
			const double coords[3] );
		      
} // end namespace PointQuery

} // end namespace FOOD

#endif // end FOOD_POINTQUERY_HPP

//---------------------------------------------------------------------------//
// end PointQuery.hpp
//---------------------------------------------------------------------------//
