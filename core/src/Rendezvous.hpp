//---------------------------------------------------------------------------//
// \file Rendezvous.hpp
// \author Stuart Slattery
// \brief Rendezvous class declaration.
//---------------------------------------------------------------------------//

#ifndef FOOD_RENDEZVOUS_HPP
#define FOOD_RENDEZVOUS_HPP

#include "TensorField.hpp"

#include <Teuchos_RCP.hpp>

#include <Tpetra_Export.hpp>

namespace FOOD
{

template<class ScalarType_T>
class Rendezvous
{
    
  public:

    //@{
    //! Typedefs.
    typedef ScalarType_T                                  ScalarType;
    typedef int                                           OrdinalType;
    typedef TensorField<ScalarType>                       TensorField_t;
    typedef Teuchos::RCP<Tensor_Field_t>                  RCP_TensorField;
    typedef Tpetra::Export<OrdinalType>                   Export_t;
    typedef Teuchos::RCP<Export_t>                        RCP_Export;
    typedef int                                           ErrorCode;
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
    RCP_Export d_domain_exporter;

    // Exporter for transfer between primary and secondary range
    // decompositions. 
    RCP_Export d_range_exporter;

  public:

    // Constructor.
    Rendezvous( RCP_TensorField domain,
		RCP_TensorField range );

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
};

} // end namespace FOOD

#endif // end FOOD_RENDEZVOUS_HPP

//---------------------------------------------------------------------------//
// end Rendezvous.hpp
//---------------------------------------------------------------------------//
