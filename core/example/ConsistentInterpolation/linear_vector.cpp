//---------------------------------------------------------------------------//
// \file cxx_main.cpp
// \author Stuart Slattery
// \brief FEMInterpolation Example 1.
//---------------------------------------------------------------------------//

#include <cassert>
#include <string>
#include <iostream>
#include <map>

#include <Domain.hpp>
#include <TensorTemplate.hpp>
#include <DFuncKernel.hpp>
#include <TensorField.hpp>
#include <FEMInterpolate.hpp>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Intrepid_FieldContainer.hpp>

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

// This example loads a fine tetrahedron mesh tagged with a 3-vector and
// interpolates it onto a coarse hexahedron mesh.
int main(int argc, char* argv[])
{
    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    if ( getDefaultComm<int>()->getRank() == 0 )
    {

    typedef Intrepid::FieldContainer<double> MDArray;

    int error;

    // The tensor template can be shared by both the range and
    // domain. 
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(1, 3, FOOD::FOOD_REAL, Teuchos::null) );

    // Need another one for the gradient. Here the gradient is rank 2 tensor.
    Teuchos::RCP<FOOD::TensorTemplate> grad_tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(2, 9, FOOD::FOOD_REAL, Teuchos::null) );

    // Set up the func_dmn mesh.
    iMesh_Instance domain_mesh;
    iMesh_newMesh("", &domain_mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    iBase_EntitySetHandle func_dmn_root_set;
    iMesh_getRootSet( domain_mesh, &func_dmn_root_set, &error );
    assert( iBase_SUCCESS == error );

    std::string domain_mesh_filename = "tagged_big_linear_tet.vtk";
    iMesh_load( domain_mesh, 
		func_dmn_root_set, 
		&domain_mesh_filename[0], 
		"", 
		&error,
		(int) domain_mesh_filename.size(),
		0 );
    assert( iBase_SUCCESS == error );

    // Set up the func_dmn mesh field.
    Teuchos::RCP<FOOD::Domain> func_dmn_domain = Teuchos::rcp(
	new FOOD::Domain(domain_mesh, func_dmn_root_set) );

    Teuchos::RCP< FOOD::DFuncKernel<double> > func_dmn_dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_TETRAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     1 ) );

    Teuchos::RCP< FOOD::TensorField<double> > func_dmn_field = Teuchos::rcp(
	new FOOD::TensorField<double>( getDefaultComm<int>(),
				       func_dmn_domain,
				       func_dmn_dfunckernel,
				       FOOD::FOOD_CARTESIAN, 
				       tensor_template,
				       Teuchos::null,
				       "FUNC_DMN_FIELD" ) );

    std::string func_dmn_tag_name = "domain";
    iBase_TagHandle func_dmn_tag;
    iMesh_getTagHandle( func_dmn_domain->getMesh(),
			&func_dmn_tag_name[0],
			&func_dmn_tag,
			&error,
			(int) func_dmn_tag_name.size() );
    assert( iBase_SUCCESS == error );

    func_dmn_field->attachToTagData( func_dmn_tag, error );
    assert( iBase_SUCCESS == error );

    // Set up the func_rng mesh.
    iMesh_Instance func_rng_mesh;
    iMesh_newMesh("", &func_rng_mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    iBase_EntitySetHandle func_rng_root_set;
    iMesh_getRootSet( func_rng_mesh, &func_rng_root_set, &error );
    assert( iBase_SUCCESS == error );

    std::string func_rng_mesh_filename = "tagged_small_99_linear_hex.vtk";
    iMesh_load( func_rng_mesh, 
		func_rng_root_set, 
		&func_rng_mesh_filename[0], 
		"", 
		&error,
		(int) func_rng_mesh_filename.size(),
		0 );
    assert( iBase_SUCCESS == error );

    // Set up the func_rng mesh field for function values.
    Teuchos::RCP<FOOD::Domain> func_rng_domain = Teuchos::rcp(
	new FOOD::Domain(func_rng_mesh, func_rng_root_set) );

    Teuchos::RCP< FOOD::DFuncKernel<double> > func_rng_dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_HEXAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     1 ) );

    Teuchos::RCP< FOOD::TensorField<double> > func_rng_field = Teuchos::rcp(
	new FOOD::TensorField<double>( getDefaultComm<int>(),
				       func_rng_domain,
				       func_rng_dfunckernel,
				       FOOD::FOOD_CARTESIAN, 
				       tensor_template,
				       Teuchos::null,
				       "FUNC_RNG_FIELD" ) );

    std::string func_rng_tag_name = "range";
    iBase_TagHandle func_rng_tag;
    iMesh_getTagHandle( func_rng_domain->getMesh(),
			&func_rng_tag_name[0],
			&func_rng_tag,
			&error,
			(int) func_rng_tag_name.size() );
    assert( iBase_SUCCESS == error );

    func_rng_field->attachToTagData( func_rng_tag, error );
    assert( iBase_SUCCESS == error );

    // Setup the gradient field.
    Teuchos::RCP< FOOD::TensorField<double> > func_rng_grad_field = Teuchos::rcp(
	new FOOD::TensorField<double>( getDefaultComm<int>(),
				       func_rng_domain,
				       func_rng_dfunckernel,
				       FOOD::FOOD_CARTESIAN, 
				       grad_tensor_template,
				       Teuchos::null,
				       "FUNC_RNG_GRAD_FIELD" ) );

    std::string func_rng_grad_tag_name = "grad_range";
    iBase_TagHandle func_rng_grad_tag;
    iMesh_getTagHandle( func_rng_domain->getMesh(),
			&func_rng_grad_tag_name[0],
			&func_rng_grad_tag,
			&error,
			(int) func_rng_grad_tag_name.size() );
    assert( iBase_SUCCESS == error );

    func_rng_grad_field->attachToTagData( func_rng_grad_tag, error );
    assert( iBase_SUCCESS == error );

    // Do interpolation.
    FOOD::FEMInterpolate<double> fem_interp( func_dmn_field, func_rng_field );
    fem_interp.setup();
    fem_interp.interpolateValueDF();

    FOOD::FEMInterpolate<double> fem_interp_grad( func_dmn_field, func_rng_grad_field );
    fem_interp_grad.setup();
    fem_interp_grad.interpolateGradDF();

    // Write the interpolated mesh to file.
    std::string interp_file = "linear_vector_output.vtk";
    iMesh_save( func_rng_domain->getMesh(),
		func_rng_domain->getMeshSet(),
		&interp_file[0],
		"",
		&error,
		(int) interp_file.size(),
		0 );
    assert( iBase_SUCCESS == error );

    } // end rank 0

    return 0;
}
