//---------------------------------------------------------------------------//
/*!
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3.0 of the License, or (at your option) any later version.
 *
 * \file cxx_main.cpp
 * \author Stuart Slattery
 * \brief Consistent interpolation example 1.
 */
//---------------------------------------------------------------------------//

#include <cassert>
#include <string>
#include <iostream>
#include <map>

#include <Domain.hpp>
#include <TensorTemplate.hpp>
#include <DFuncKernel.hpp>
#include <TensorField.hpp>
#include <ConsistentScheme.hpp>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ENull.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_CommHelpers.hpp>

template<class Ordinal>
Teuchos::RCP<const Teuchos::Comm<Ordinal> > getDefaultComm()
{
#ifdef HAVE_MPI
    return Teuchos::DefaultComm<Ordinal>::getComm();
#else
    return Teuchos::rcp(new Teuchos::SerialComm<Ordinal>() );
#endif
}

// This consistent interpolation example loads a linear quad mesh
// tagged with a scalar field and interpolates it along with the gradient
// pullback onto a delaunay tri mesh mesh. (Function domain = func_dmn,
// function range = func_rng ) 
int main(int argc, char* argv[])
{
    // Setup communication.
    Teuchos::GlobalMPISession mpiSession(&argc,&argv);
    Teuchos::RCP<const Teuchos::Comm<int> > comm = 
	Teuchos::DefaultComm<int>::getComm();

    if ( getDefaultComm<int>()->getRank() == 0 ) // Force scalar execution.
    {

    int error;

    // The tensor template can be shared by both the range and
    // domain value. 
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL, Teuchos::null) );

    // Need another one for the gradient. Here the gradient is a 3-vector.
    Teuchos::RCP<FOOD::TensorTemplate> grad_tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(1, 3, FOOD::FOOD_REAL, Teuchos::null) );

    // Set up the func_dmn mesh.
    iMesh_Instance func_dmn_mesh;
    iMesh_newMesh("", &func_dmn_mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    iBase_EntitySetHandle func_dmn_root_set;
    iMesh_getRootSet( func_dmn_mesh, &func_dmn_root_set, &error );
    assert( iBase_SUCCESS == error );

    std::string func_dmn_mesh_filename = "tagged_quad_flat_surf.vtk";
    iMesh_load( func_dmn_mesh, 
		func_dmn_root_set, 
		&func_dmn_mesh_filename[0], 
		"", 
		&error,
		(int) func_dmn_mesh_filename.size(),
		0 );
    assert( iBase_SUCCESS == error );

    // Set up the func_dmn mesh field.
    Teuchos::RCP<FOOD::Domain> func_dmn_domain = Teuchos::rcp(
	new FOOD::Domain(func_dmn_mesh, func_dmn_root_set, FOOD::FOOD_MBCN) );

    Teuchos::RCP< FOOD::DFuncKernel<double> > func_dmn_dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_FACE,
						     iMesh_QUADRILATERAL,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
						     2 ) );

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

    std::string func_rng_mesh_filename = "tagged_delaunay_flat_surf.vtk";
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
	new FOOD::Domain(func_rng_mesh, func_rng_root_set, FOOD::FOOD_MBCN) );

    Teuchos::RCP< FOOD::DFuncKernel<double> > func_rng_dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_FACE,
						     iMesh_TRIANGLE,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     FOOD::FOOD_SHARDSCN,
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
    FOOD::ConsistentScheme<double> fem_interp_val( func_dmn_field, 
						   func_rng_field );
    fem_interp_val.setup();
    fem_interp_val.transferValueDF();

    FOOD::ConsistentScheme<double> fem_interp_grad( func_dmn_field, 
						    func_rng_grad_field );
    fem_interp_grad.setup();
    fem_interp_grad.transferGradDF();

    // Write the interpolated mesh to file.
    std::string interp_file = "flat_surf_output.vtk";
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
