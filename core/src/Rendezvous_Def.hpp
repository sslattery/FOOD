//---------------------------------------------------------------------------//
// \file Rendezvous_Def.hpp
// \author Stuart Slattery
// \brief Rendezvous class definition.
//---------------------------------------------------------------------------//

#ifndef FOOD_RENDEZVOUS_DEF_HPP
#define FOOD_RENDEZVOUS_DEF_HPP

#include <cassert>

#include <Tpetra_Vector.hpp>

namespace FOOD
{
 
/*!
 * \brief Constructor.
 */
template<class ScalarType>
Rendezvous<ScalarType>::Rendezvous(RCP_TensorField domain, 
				   RCP_TensorField range )
    : d_domain_primary(domain)
    , d_range_primary(range)
    , d_domain_secondary(0)
    , d_range_secondary(0)
    , d_domain_exporter(0)
    , d_range_exporter(0)
{
    assert( domain->getTensorFieldTemplate() == 
	    range->getTensorFieldTemplate() );
}

template<class ScalarType>
Rendezvous<ScalarType>::~Rendezvous()
{ /* ... */ }

} // end namespace FOOD

#endif // end FOOD_RENDEZVOUS_DEF_HPP

//---------------------------------------------------------------------------//
// end Rendezvous_Def.hpp
//---------------------------------------------------------------------------//
