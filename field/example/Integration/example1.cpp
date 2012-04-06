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
// Integrate the field applied to a mesh. The field data is of type double so
// all objects below are templated with double as the type.
int main(int argc, char* argv[])
{
    int error;

    // Create a tensor template. This describes the mathematical properties of
    // the tensor we will create. Here, the tensor (i.e. the tag data) is a
    // scalar quantity (arg1 = 0), each tag has one scalar (arg2 = 1), and
    // tensor is composed of real numbers (arg3 = FOOD::FOOD_REAL). Also note
    // that tensor_template is actually a smart pointer to a
    // TensorTemplate. We just wrap Teuchos::rcp() around the new call when
    // allocating the object on the heap.
    Teuchos::RCP<FOOD::TensorTemplate> tensor_template = Teuchos::rcp(
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL) );

    // Set up the mesh. Use iMesh here to load the mesh into the root set.
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

    // Create a domain for the field. This is the set in the mesh that the
    // field is applied to. We are operating on the mesh we just loaded from
    // file (arg1 = mesh), the set we want to apply the field to is the root
    // set (arg2 = root_set), and the mesh was created using moab based tools
    // so we expect the mesh entities to use the moab canonical numbering
    // scheme for their connectivity (arg3 = FOOD::FOOD_MBCN)
    Teuchos::RCP<FOOD::Domain> domain = Teuchos::rcp(
	new FOOD::Domain(mesh, root_set, FOOD::FOOD_MBCN) );

    // Setup a kernel factory. This will automatically build the distribution
    // function kernels for us (i.e. the basis).
    FOOD::DFuncKernelFactory<double> kernel_factory;

    // Use the kernel factory to automatically create a basis. In this
    // example, we are applying the basis to the hexahedrons in the mesh (arg1
    // = iBase_REGION, arg2 = iMesh_HEXAHEDRON), we want a finite element
    // basis (arg3 = FOOD::FOOD_FEM), we will be operating in the Hilbert
    // space spanned by the gradient operator, a good default if you only care
    // about field values and not operator values (arg4 = FOOD:FOOD_HGRAD),
    // and we want a quadratic basis as the mesh in this example consists of
    // quadratic hexahedron elements (arg5 = 2)
    Teuchos::RCP< FOOD::DFuncKernel<double> > dfunckernel = 
	kernel_factory.create( iBase_REGION,
			       iMesh_HEXAHEDRON,
			       FOOD::FOOD_FEM,
			       FOOD::FOOD_HGRAD,
			       2 );

    // Now create the actual tensor field that will apply the discretization
    // to the mesh and its data. We are applying this field to the domain we
    // specified (arg1 = domain), we want to apply the discretization we just
    // created with the kernel factory (arg2 = dfunckernel), we want this
    // tensor to have the mathematical properties of the tensor template we
    // defined (arg3 = tensor_template), and we'll call this field EX1_FIELD
    // (arg4 = "EX1_FIELD")
    Teuchos::RCP< FOOD::TensorField<double> > field = Teuchos::rcp(
	new FOOD::TensorField<double>( domain,
				       dfunckernel,
				       tensor_template,
				       "EX1_FIELD" ) );

    // Now we need to apply the discretization to the mesh data. First, we get
    // the tag for the data we want. In this example, the mesh has been tagged
    // with a tag called "domain"
    std::string tag_name = "domain";
    iBase_TagHandle tag;
    iMesh_getTagHandle( domain->getMesh(),
			&tag_name[0],
			&tag,
			&error,
			(int) tag_name.size() );
    assert( iBase_SUCCESS == error );

    // We call attachToTagData to apply the tag data to the field. Here I'm
    // returning the error code but don't expect this in the future as FOOD's
    // error handling policy is to throw C++ exceptions. I need to write an
    // ITAPS execption class that correctly throws itaps error codes as
    // exceptions. 
    field->attachToTagData( tag, error );
    assert( iBase_SUCCESS == error );

    // Call integrateCells to integrate the nodal data. The resulting integral
    // will be applied to the hexahedron volumes in a tag called
    // EX1_FIELD_INTEGRAL
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
