//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file DFuncKernelFactory_Def.hpp
 * \author Stuart Slattery
 * \brief Factory method defintion for basis generation.
 */
//---------------------------------------------------------------------------//

#ifndef FOOD_DFUNCKERNELFACTORY_DEF_HPP
#define FOOD_DFUNCKERNELFACTORY_DEF_HPP

#include <cassert>

#include "DiscretizationTypes.hpp"
#include "IntrepidKernel.hpp"

#include <iBase.h>
#include <iMesh.h>

#include <Intrepid_FieldContainer.hpp>
#include <Intrepid_Basis.hpp>
#include <Intrepid_HCURL_HEX_I1_FEM.hpp>
#include <Intrepid_HCURL_QUAD_I1_FEM.hpp>
#include <Intrepid_HCURL_TET_I1_FEM.hpp>
#include <Intrepid_HCURL_TRI_I1_FEM.hpp>
#include <Intrepid_HDIV_HEX_I1_FEM.hpp>
#include <Intrepid_HDIV_QUAD_I1_FEM.hpp>
#include <Intrepid_HDIV_TET_I1_FEM.hpp>
#include <Intrepid_HDIV_TRI_I1_FEM.hpp>
#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid_HGRAD_HEX_C2_FEM.hpp>
#include <Intrepid_HGRAD_LINE_C1_FEM.hpp>
#include <Intrepid_HGRAD_QUAD_C1_FEM.hpp>
#include <Intrepid_HGRAD_QUAD_C2_FEM.hpp>
#include <Intrepid_HGRAD_TET_C1_FEM.hpp>
#include <Intrepid_HGRAD_TET_C2_FEM.hpp>
#include <Intrepid_HGRAD_TRI_C1_FEM.hpp>
#include <Intrepid_HGRAD_TRI_C2_FEM.hpp>
#include <Intrepid_HGRAD_WEDGE_C1_FEM.hpp>
#include <Intrepid_HGRAD_WEDGE_C2_FEM.hpp>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar>
DFuncKernelFactory<Scalar>::DFuncKernelFactory()
{ /* ... */ }

/*!
 * \brief Destructor.
 */
template<class Scalar>
DFuncKernelFactory<Scalar>::~DFuncKernelFactory()
{ /* ... */ }

/*!
 * \brief Factory Method
 */
template<class Scalar>
Teuchos::RCP< DFuncKernel<Scalar> > 
DFuncKernelFactory<Scalar>::create( const int entity_topology,
				    const int discretization_type,
				    const int function_space_type,
				    const int degree )

{
    typedef typename IntrepidKernel<Scalar>::MDArray MDArray;
    Teuchos::RCP< DFuncKernel<Scalar> > new_dfunckernel;
    
    switch ( function_space_type )
    {
	
	case FOOD_HGRAD:

	    switch ( entity_topology )
	    {
		case iMesh_LINE_SEGMENT:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp( 
					    new Intrepid::Basis_HGRAD_LINE_C1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
				    );
			    }
			    else
			    {
				assert( degree == 1 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		case iMesh_TRIANGLE:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HGRAD_TRI_C1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
				    );
			    }
			    else if ( degree == 2 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HGRAD_TRI_C2_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
				    );
			    }
			    else
			    {
				assert( degree == 1 ||
					degree == 2 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;

		case iMesh_QUADRILATERAL:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
				    );
			    }
			    else if ( degree == 2 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 ||
					degree == 2 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;

		case iMesh_TETRAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>( 
					Teuchos::rcp(
					    new Intrepid::Basis_HGRAD_TET_C1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else if ( degree == 2 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HGRAD_TET_C2_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 ||
					degree == 2 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;
		    
		case iMesh_HEXAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HGRAD_HEX_C1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else if ( degree == 2 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HGRAD_HEX_C2_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 ||
					degree == 2 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;

		default:

		    assert( iMesh_LINE_SEGMENT  == entity_topology ||
			    iMesh_TRIANGLE      == entity_topology ||
			    iMesh_QUADRILATERAL == entity_topology ||
			    iMesh_TETRAHEDRON   == entity_topology ||
			    iMesh_HEXAHEDRON    == entity_topology );
	    }
	
	    break;

	case FOOD_HDIV:

	    switch ( entity_topology )
	    {

		case iMesh_TRIANGLE:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HDIV_TRI_I1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
				    );
			    }
			    else
			    {
				assert( degree == 1 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;

		case iMesh_QUADRILATERAL:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HDIV_QUAD_I1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;

		case iMesh_TETRAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HDIV_TET_I1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;
		    
		case iMesh_HEXAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HDIV_HEX_I1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;

		default:

		    assert( iMesh_TRIANGLE      == entity_topology ||
			    iMesh_QUADRILATERAL == entity_topology ||
			    iMesh_TETRAHEDRON   == entity_topology ||
			    iMesh_HEXAHEDRON    == entity_topology );
	    }

	    break;

	case FOOD_HCURL:

	    switch ( entity_topology )
	    {
		case iMesh_TRIANGLE:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HCURL_TRI_I1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;

		case iMesh_QUADRILATERAL:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HCURL_QUAD_I1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;

		case iMesh_TETRAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HCURL_TET_I1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;
		    
		case iMesh_HEXAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( degree == 1 )
			    {
				new_dfunckernel = Teuchos::rcp( 
				    new IntrepidKernel<Scalar>(
					Teuchos::rcp(
					    new Intrepid::Basis_HCURL_HEX_I1_FEM<Scalar,MDArray>() ),
					    entity_topology,
					    discretization_type,
					    function_space_type )
					);
			    }
			    else
			    {
				assert( degree == 1 );
			    }

			    break;

			default:

			    assert ( FOOD_FEM == discretization_type );
		    }

		    break;

		default:

		    assert( iMesh_TRIANGLE      == entity_topology ||
			    iMesh_QUADRILATERAL == entity_topology ||
			    iMesh_TETRAHEDRON   == entity_topology ||
			    iMesh_HEXAHEDRON    == entity_topology );
	    }

	    break;

	default:

	    assert( FOOD_HGRAD == function_space_type ||
		    FOOD_HDIV  == function_space_type ||
		    FOOD_HCURL == function_space_type );
    }

    return new_dfunckernel;
}

} // end namepsace FOOD
	
#endif // end FOOD_DFUNCKERNELFACTORY_DEF_HPP

//---------------------------------------------------------------------------//
// end DFuncKernelFactory_Def.hpp
//---------------------------------------------------------------------------//
