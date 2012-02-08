//---------------------------------------------------------------------------//
// \file cxx_main.cpp
// \author Stuart Slattery
// \brief Driver for tet part interpolation example.
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
	new FOOD::TensorTemplate(0, 1, FOOD::FOOD_REAL, Teuchos::null) );

    // Set up the fine mesh.
    iMesh_Instance fine_mesh;
    iMesh_newMesh("", &fine_mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    iBase_EntitySetHandle fine_root_set;
    iMesh_getRootSet( fine_mesh, &fine_root_set, &error );
    assert( iBase_SUCCESS == error );

    std::string fine_mesh_filename = "tagged_fine.vtk";
    iMesh_load( fine_mesh, 
		fine_root_set, 
		&fine_mesh_filename[0], 
		"", 
		&error,
		(int) fine_mesh_filename.size(),
		0 );
    assert( iBase_SUCCESS == error );

    // Set up the fine mesh field.
    Teuchos::RCP<FOOD::Domain> fine_domain = Teuchos::rcp(
	new FOOD::Domain(fine_mesh, fine_root_set) );

    Teuchos::RCP< FOOD::DFuncKernel<double> > fine_dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_TETRAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     1 ) );

    Teuchos::RCP< FOOD::TensorField<double> > fine_field = Teuchos::rcp(
	new FOOD::TensorField<double>( getDefaultComm<int>(),
				       fine_domain,
				       fine_dfunckernel,
				       FOOD::FOOD_CARTESIAN, 
				       tensor_template,
				       Teuchos::null,
				       "FINE_FIELD" ) );

    std::string fine_tag_name = "domain";
    iBase_TagHandle fine_tag;
    iMesh_getTagHandle( fine_domain->getMesh(),
			&fine_tag_name[0],
			&fine_tag,
			&error,
			(int) fine_tag_name.size() );
    assert( iBase_SUCCESS == error );

    fine_field->attachToTagData( fine_tag, error );
    assert( iBase_SUCCESS == error );

    // Set up the coarse mesh.
    iMesh_Instance coarse_mesh;
    iMesh_newMesh("", &coarse_mesh, &error, 0);
    assert( iBase_SUCCESS == error );

    iBase_EntitySetHandle coarse_root_set;
    iMesh_getRootSet( coarse_mesh, &coarse_root_set, &error );
    assert( iBase_SUCCESS == error );

    std::string coarse_mesh_filename = "tagged_coarse_99_hex.vtk";
    iMesh_load( coarse_mesh, 
		coarse_root_set, 
		&coarse_mesh_filename[0], 
		"", 
		&error,
		(int) coarse_mesh_filename.size(),
		0 );
    assert( iBase_SUCCESS == error );

    // Set up the coarse mesh field for function values.
    Teuchos::RCP<FOOD::Domain> coarse_domain = Teuchos::rcp(
	new FOOD::Domain(coarse_mesh, coarse_root_set) );

    Teuchos::RCP< FOOD::DFuncKernel<double> > coarse_dfunckernel =
	Teuchos::rcp( new FOOD::DFuncKernel<double>( iBase_REGION,
						     iMesh_HEXAHEDRON,
						     iBase_VERTEX,
						     iMesh_POINT,
						     FOOD::FOOD_CARTESIAN,
						     FOOD::FOOD_FEM,
						     FOOD::FOOD_HGRAD,
						     1 ) );

    Teuchos::RCP< FOOD::TensorField<double> > coarse_field = Teuchos::rcp(
	new FOOD::TensorField<double>( getDefaultComm<int>(),
				       coarse_domain,
				       coarse_dfunckernel,
				       FOOD::FOOD_CARTESIAN, 
				       tensor_template,
				       Teuchos::null,
				       "COARSE_FIELD" ) );

    std::string coarse_tag_name = "range";
    iBase_TagHandle coarse_tag;
    iMesh_getTagHandle( coarse_domain->getMesh(),
			&coarse_tag_name[0],
			&coarse_tag,
			&error,
			(int) coarse_tag_name.size() );
    assert( iBase_SUCCESS == error );

    coarse_field->attachToTagData( coarse_tag, error );
    assert( iBase_SUCCESS == error );

    // Do interpolation.
    FOOD::FEMInterpolate<double> fem_interp( fine_field, coarse_field );
    fem_interp.setup();
    fem_interp.interpolateValueDF();

    // Write the interpolated mesh to file.
    std::string interp_file = "interpolated_coarse_99.vtk";
    iMesh_save( coarse_domain->getMesh(),
		coarse_domain->getMeshSet(),
		&interp_file[0],
		"",
		&error,
		(int) interp_file.size(),
		0 );
    assert( iBase_SUCCESS == error );

    } // end rank 0

    return 0;
}
