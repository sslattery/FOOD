//---------------------------------------------------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// \file Rendezvous.hpp
// \author Stuart Slattery
// \brief Rendezvous class declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_RENDEZVOUS_HPP
#define FOOD_RENDEZVOUS_HPP

#include "TensorField.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>

#include <Tpetra_Export.hpp>

namespace FOOD
{

template<class Scalar>
class Rendezvous
{
    
  public:

    //@{
    //! Typedefs.
    typedef int                                           OrdinalType;
    typedef TensorField<Scalar>                           TensorField_t;
    typedef Teuchos::RCP<TensorField_t>                   RCP_TensorField;
    typedef Tpetra::Export<OrdinalType>                   Export_t;
    typedef Teuchos::RCP<Export_t>                        RCP_Export;
    typedef int                                           ErrorCode;
    typedef iBase_EntityHandle                            EntityHandle;
    //@}

  private:

    // Domain with primary decomposition.
    RCP_TensorField d_domain_primary;

    // Range with primary decomposition.
    RCP_TensorField d_range_primary;

    // Domain with secondary decomposition.
    RCP_TensorField d_domain_secondary;

    // Range with secondary decomposition.
    RCP_TensorField d_range_secondary;

    // Exporter for transfer between primary and secondary domain
    // decompositions. 
    RCP_Export d_domain_export;

    // Exporter for transfer between primary and secondary range
    // decompositions. 
    RCP_Export d_range_export;

  public:

    // Constructor.
    Rendezvous( RCP_TensorField domain,	RCP_TensorField range );

    // Destructor.
    ~Rendezvous();

    // Do parallel rendezvous to generate secondary decompositions.
    void createSecondaryDecompositions();

    // Copy the domain degrees of freedom from the primary decomposition to
    // the secondary. 
    void domainCopyPrimaryToSecondary();

    // Copy the range degrees of freedom from the secondary decomposition to
    // the primary.
    void rangeCopySecondaryToPrimary();

    //! Get the domain primary decomposition.
    RCP_TensorField getRendezvousDomainPrimary() const
    { return d_domain_primary; }

    //! Get the range primary decomposition.
    RCP_TensorField getRendezvousRangePrimary() const
    { return d_range_primary; }

    //! Get the domain secondary decomposition.
    RCP_TensorField getRendezvousDomainSecondary() const
    { return d_domain_secondary; }

    //! Get the range secondary decomposition.
    RCP_TensorField getRendezvousRangeSecondary() const
    { return d_domain_secondary; }

  private:

    // Compute a global axis aligned bounding box that bounds the intersection
    // of the domain and range mesh sets. 
    Teuchos::Tuple<double,6> computeIntersectionBoundingBox();

    // Intersection test for two axis aligned boxes. Return the resulting
    // intersecting axis aligned box.
    Teuchos::Tuple<double,6> boxIntersectionTest( Teuchos::Tuple<double,6> box_A,
						  Teuchos::Tuple<double,6> box_B );
};

} // end namespace FOOD

#include "Rendezvous_Def.hpp"

#endif // end FOOD_RENDEZVOUS_HPP

//---------------------------------------------------------------------------//
// end Rendezvous.hpp
//---------------------------------------------------------------------------//
