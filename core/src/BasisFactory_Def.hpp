//---------------------------------------------------------------------------//
// \file BasisFactory_Def.hpp
// \author Stuart Slattery
// \brief Factory method defintion for basis generation.
//---------------------------------------------------------------------------//

#ifndef FOOD_BASISFACTORY_DEF_HPP
#define FOOD_BASISFACTORY_DEF_HPP

#include <cassert>

namespace FOOD
{

/*!
 * \brief Constructor.
 */
template<class Scalar, class ArrayScalar>
BasisFactory<Scalar,ArrayScalar>::BasisFactory()
{ /* ... */ }

/*!
 * \brief Destructor.
 */
template<class Scalar, class ArrayScalar>
BasisFactory<Scalar,ArrayScalar>::~BasisFactory()
{ /* ... */ }

/*!
 * \brief Factory Method
 */
template<class Scalar, class ArrayScalar>
Teuchos::RCP< Intrepid::Basis<Scalar,ArrayScalar> > 
BasisFactory<Scalar,ArrayScalar>::create( const int entity_topology,
					  const int discretization_type,
					  const int basis_function_space,
					  const int basis_degree )

{
    Teuchos::RCP< Intrepid::Basis<Scalar,ArrayScalar> > new_basis;
    
    switch ( basis_function_space )
    {
	
	case FOOD_HGRAD:

	    switch ( entity_topology )
	    {
		case iMesh_LINE_SEGMENT:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HGRAD_LINE_C1_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		case iMesh_TRIANGLE:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HGRAD_TRI_C1_FEM<Scalar,ArrayScalar>() );
			    }
			    else if ( basis_degree == 2 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HGRAD_TRI_C2_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;

		case iMesh_QUADRILATERAL:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Scalar,ArrayScalar>() );
			    }
			    else if ( basis_degree == 2 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;

		case iMesh_TETRAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HGRAD_TET_C1_FEM<Scalar,ArrayScalar>() );
			    }
			    else if ( basis_degree == 2 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HGRAD_TET_C2_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;
		    
		case iMesh_HEXAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HGRAD_HEX_C1_FEM<Scalar,ArrayScalar>() );
			    }
			    else if ( basis_degree == 2 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HGRAD_HEX_C2_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;
	    }
	
	    break;

	case FOOD_HDIV:

	    switch ( entity_topology )
	    {

		case iMesh_TRIANGLE:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HDIV_TRI_I1_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;

		case iMesh_QUADRILATERAL:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HDIV_QUAD_I1_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;

		case iMesh_TETRAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HDIV_TET_I1_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;
		    
		case iMesh_HEXAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HDIV_HEX_I1_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;
	    }

	    break;

	case FOOD_HCURL:

	    switch ( entity_topology )
	    {
		case iMesh_TRIANGLE:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HCURL_TRI_I1_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;

		case iMesh_QUADRILATERAL:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HCURL_QUAD_I1_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;

		case iMesh_TETRAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HCURL_TET_I1_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;
		    
		case iMesh_HEXAHEDRON:

		    switch ( discretization_type )
		    {
			case FOOD_FEM:

			    if ( basis_degree == 1 )
			    {
				new_basis = Teuchos::rcp( 
				    new Intrepid::Basis_HCURL_HEX_I1_FEM<Scalar,ArrayScalar>() );
			    }

			    break;
		    }

		    break;
	    }

	    break;
    }

    return new_basis;
}

} // end namepsace FOOD
	
#endif // end FOOD_BASISFACTORY_DEF_HPP

//---------------------------------------------------------------------------//
// end BasisFactory_Def.hpp
//---------------------------------------------------------------------------//
