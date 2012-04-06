//---------------------------------------------------------------------------//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// \file cxx_main.cpp
// \author Stuart Slattery
// \brief Integration example 1.
//---------------------------------------------------------------------------//

#include <cassert>
#include <string>
#include <iostream>
#include <map>

#include <Domain.hpp>
#include <TensorTemplate.hpp>
#include <DFuncKernel.hpp>
#include <DFuncKernelFactory.hpp>
#include <TensorField.hpp>

#include <iBase.h>
#include <iMesh.h>

#include <Teuchos_RCP.hpp>

//---------------------------------------------------------------------------//
// Integrate the field applied to a mesh.
int main(int argc, char* argv[])
{
    int error;

    // Setup a kernel factory.
    FOOD::DFuncKernelFactory<double> kernel_factory;

    // The tensor template can be shared by both the range and
    // domain. 
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL) );

    // Set up the mesh.
    iMesh_Instance mesh;
    iMesh_newMesh("", &mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    iBase_EntitySetHandle root_set;
    iMesh_getRootSet( mesh, &root_set, &error );
    assert( iBase_SUCCESS == error );

    std::string mesh_filename = "hex_domain.vtk";
    iMesh_load( mesh, 
		root_set, 
		&mesh_filename[0], 
		"", 
		&error,
		(int) mesh_filename.size(),
		0 );
    assert( iBase_SUCCESS == error );

    // Set up the field.
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel = 
	kernel_factory.create( iBase_REGION,
			       iMesh_HEXAHEDRON,
			       FOOD::FOOD_FEM,
			       FOOD::FOOD_HGRAD,
			       2 );

    Teuchos::RCP< FOOD::TensorField<double> > field = Teuchos::rcp(
	new FOOD::TensorField<double>( domain,
				       dfunckernel,
				       tensor_template,
				       "EX1_FIELD" ) );

    std::string tag_name = "domain";
    iBase_TagHandle tag;
    iMesh_getTagHandle( domain->getMesh(),
			&tag_name[0],
			&tag,
			&error,
			(int) tag_name.size() );
    assert( iBase_SUCCESS == error );

    field->attachToTagData( tag, error );
    assert( iBase_SUCCESS == error );

    // Integrate the degrees of freedom over the cells.
    field->integrateCells();

    // Write the mesh to file.
    std::string out_file = "example1_output.vtk";
    iMesh_save( domain->getMesh(),
		domain->getMeshSet(),
		&out_file[0],
		"",
		&error,
		(int) out_file.size(),
		0 );
    assert( iBase_SUCCESS == error );

    return 0;
}
