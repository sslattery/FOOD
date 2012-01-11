#ifndef iField_H
#define iField_H

#include "iBase.h"
#include "iMesh.h"
#include "iGeom.h"
#include "iRel.h"
    
#ifdef __cplusplus
extern "C" {
#endif

enum iField_Precision {
    iField_Precision_MIN = 0,
        /**< MIN symbol used to facilitate iteration over enum values */
    iField_BOOLEAN = iField_Precision_MIN,
        /**< Description unknown. */
    iField_UCHAR,
        /**< Description unknown. */
    iField_INTEGER,
        /**< Description unknown. */
    iField_FLOAT,
        /**< Description unknown. */
    iField_DOUBLE,
        /**< Description unknown. */
    iField_QUAD,
        /**< Description unknown. */
    iField_Precision_MAX = iField_QUAD
        /**< MAX symbol used to facilitate iteration over enum values */
};

enum iField_AlgType {
    iField_AlgType_MIN = 0,
        /**< MIN symbol used to facilitate iteration over enum values */
    iField_LOGICAL = iField_AlgType_MIN,
        /**< Description unknown. */
    iField_INTEGRAL,
        /**< Description unknown. */
    iField_REAL,
        /**< Description unknown. */
    iField_COMPLEX,
        /**< Description unknown. */
    iField_AlgType_MAX = iField_COMPLEX
        /**< MAX symbol used to facilitate iteration over enum values */
};

enum iField_CoordType {
    iField_CoordType_MIN = 0,
        /**< MIN symbol used to facilitate iteration over enum values */
    iField_CARTESIAN = iField_CoordType_MIN,
        /**< Description unknown. */
    iField_CYLINDRICAL,
        /**< Description unknown. */
    iField_SPHERICAL,
        /**< Description unknown. */
    iField_CoordType_MAX = iField_SPHERICAL
        /**< MAX symbol used to facilitate iteration over enum values */
};

enum iField_StorageHint {
    iField_StorageHint_MIN = 0,
        /**< MIN symbol used to facilitate iteration over enum values */
    iField_BLOCKED = iField_StorageHint_MIN,
        /**< Description unknown. */
    iField_INTERLEAVED,
        /**< Description unknown. */
    iField_MIXED,
        /**< Description unknown. */
    iField_PER_ENTITY,
        /**< Description unknown. */
    iField_StorageHint_MAX = iField_PER_ENTITY
        /**< MAX symbol used to facilitate iteration over enum values */
};
  
typedef void* iField_Instance;
typedef struct iField_DomainHandle_Private* iField_DomainHandle;
typedef struct iField_TensorHandle_Private* iField_TensorHandle;
typedef struct iField_TensorTemplateHandle_Private* iField_TensorTemplateHandle;
typedef struct iField_QuantHandle_Private* iField_QuantHandle;
typedef struct iField_UnitsHandle_Private* iField_UnitsHandle;
typedef struct iField_Handle_Private* iField_Handle;

/***************************************************************************//**
 * \ingroup iField_Initialization
 * \brief Create an iField instance.
 *
 * In typical usage, the fields associated with a single iMesh instance
 * (and its various subsets) will be part of the same Field instance.
 * Question: do we want to enforce this typical usage by requiring an
 * iMesh instance as an argument in this constructor?  If so, do we also
 * want to pass an iRel instance so that that connection will be
 * mediated "properly" through iRel?  Or should we leave that up to the
 * application to handle properly?
 ******************************************************************************/
void iField_create(
    const char* options,
        /**< [in] Implementation-specific options. */
    iField_Instance* instance,
        /**< [in] iField instance handle */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    const int options_len
        /**< [in] Length of the option string */
);

/***************************************************************************//**
 * \ingroup iField_Initialization
 * \brief Destroy a field interface object
 *
 * Naturally, this call also destroys all field data stored in the
 * interface object.
 ******************************************************************************/
void iField_destroy(
    iField_Instance instance,
        /**< [in] iField instance handle */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_ErrorHandling
 * \brief Get a description of the error returned from the last iField function
 *
 * Get a description of the error returned from the last iField
 * function. Space must be allocated by the caller.
 ******************************************************************************/
void iField_getDescription(
    iField_Instance instance,
        /**< [in] iField instance handle */
    char* descr,
        /**< [in,out] Pointer to a character string to be filled with a
        description of the error from the last iField function */
    int descr_len
        /**< [in] Length of the character string pointed to by descr */
);

/***************************************************************************//**
 * \ingroup iField_TensorTemplates
 * \brief Create a tensor template; return its handle.
 *
 * Create a tensor template.  Once a tensor template is set up, its
 * basic data can be read, but not written.  So, for example, if you
 * want to change from a 2-vector to a 3-vector, you need to have
 * another tensor template for that tensor.
 ******************************************************************************/
void iField_createTensorTemplate(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    int order,
        /**< Tensor order */
    int num_comp,
        /**< Number of components */
    int alg_type,
        /**< [in] Is this data real, complex, etc? (enum) */
    iField_QuantHandle quant_handle,
        /**< [in] Physical quantity this tensor represents. */
    iField_TensorTemplateHandle* handle,
        /**< [out] The newly created tensor template. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_TensorTemplates
 * \brief Retrieve the order of this tensor template.
 *
 *
 ******************************************************************************/
void iField_getTensorTemplateOrder(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorTemplateHandle handle,
        /**< [in] The tensor template being queried. */
    int* order,
        /**< [out] The order of this tensor template. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_TensorTemplates
 * \brief Retrieve the number of components for this tensor template.
 *
 *
 ******************************************************************************/
void iField_getTensorTemplateNumComponents(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorTemplateHandle handle,
        /**< [in] The tensor template being queried. */
    int* num_comp,
        /**< [out] # of components of this tensor template. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_TensorTemplates
 * \brief Retrieve the algebraic type of this tensor template.
 *
 *
 ******************************************************************************/
void iField_getTensorTemplateAlgType(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorTemplateHandle handle,
        /**< [in] The tensor template being queried. */
    int* alg_type,
        /**< [out] The algebraic type of this tensor template. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_TensorTemplates
 * \brief Retrieve the physical quantity this tensor template represents.
 *
 ******************************************************************************/
void iField_getTensorTemplateQuantity(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorTemplateHandle handle,
        /**< [in] The tensor template being queried. */
    iField_QuantHandle* quant_handle,
        /**< [out] The quantity of this tensor template. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_TensorTemplates
 * \brief Destroy a tensor template.
 *
 * This should refuse to destroy if the template is in use.
 ******************************************************************************/
void iField_destroyTensorTemplate(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorTemplateHandle handle,
        /**< [in] Tensor template to destroy. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_TensorMembers
 * \brief Create a tensor field from a template, distribution
 *  functions, units, and domain.
 *
 ******************************************************************************/
void iField_createTensorField(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_DomainHandle domain_handle,
        /**< [in] The domain this field is defined on. */
    int precision,
        /**< [in] Is this data single precision, double precision, etc? (enum) */
    iField_DFuncKernel dfunc_kernel,
        /**< [in] A function callback for converting spatial coordinates into
        the value for a distribution function. */
    iField_TensorTemplateHandle tt_handle,
        /**< [in] The tensor template specifies whether this is, for instance,
        a velocity, or a temperature, or a stress. */
    iField_UnitsHandle units_handle,
        /**< [in] SI, cgs, slugs-ft-sec, or user-defined. */
    const char* name,
        /**< [in] A string giving a name for this tensor. */
    iField_TensorHandle* tensor_handle,
        /**< [out] The new tensor field. */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    const int name_len
        /**< [in] Length of name string. */
);

/***************************************************************************//**
 * \ingroup iField_TensorMembers
 * \brief Create multiple tensor fields from a template, distribution
 *
 *  functions, units, and coordinate system.
 *  The concept this multi-field create method provides are...
 *    1. hint to implementation of storing them together,
 *    2. possible convenience for caller to set all dofs together
 *    3. storage is predictable/computable because all fields obey
 *    properties that make this possible. This also has implications in
 *    array-based methods of getting/setting dofs as well.
 *  Have to look into Fortran size requirements here!
 * Restrictions:
 *  1. The resulting component fields can be evaluated and
 *  set/getDOF'ed individually, but can't be destroyed individually.
 *  2. All fields must use the same dfunc and data storage precision.
 *  3. Units need not be consistent between fields, as they're each
 *  defined separately.
 ******************************************************************************/
void iField_createTensorFields(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_DomainHandle domain_handle,
        /**< [in] The domain this field is defined on. */
    int precision,
        /**< [in] Is this data single precision, double precision, etc? (enum) */
    iField_DFuncKernel dfunc_kernel,
        /**< [in] A function callback for converting spatial coordinates into the
        value for a distribution function. */
    iField_TensorTemplateHandle* tt_handles,
        /**< [in] One tensor template for each tensor to be created.  */
    iField_UnitsHandle* units_handles,
        /**< [in] Units for each tensor to be created. */
    char** names,
        /**< [in] An array of strings giving names of tensors (one more name
        than templates, so that the compound will also have a name). */
    const int num_tensors,
        /**< [in] How many new tensors? */
    iField_TensorHandle* compound_handle,
        /**< [out] A separate handle for the collection of tensors.  */
    iField_TensorHandle** tensor_handles,
        /**< [in,out] The new tensor fields. */
    int* tensor_handles_allocated,
        /**< [in,out] Allocated size of tensor handle array. */
    int* tensor_handles_size,
        /**< [out] Number of tensors created and returned. */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    const int* name_len
        /**< Lengths of name strings. [CFOG */
);

/***************************************************************************//**
 * \ingroup iField_DOFStorage
 * \brief Indicate that the given field should use tag data stored in
 *  the given taghandle.
 *
 * This function explicitly assumes that the tag data will be used as
 * the field DOFs; any existing DOFs will be destroyed by this call in
 * favor of retaining the values stored in the tag.
 * iField set/get dofs will have the same effect as iMesh set/get tags.
 * A default value is provided so that any entities that should have
 * DOFs but don't can be assigned values.  This value -must- be set to
 * tags (instead of held back as a default return value from getDOFs) so
 * that getTag will properly return the field dofs.
 * A (simple or compound) field that is stored as a tag can not respond
 * to suggestions for storage.
 * All floating point field types will be converted to double precision
 * by this call, including a change in the field precision.
 * Floating point fields must attach to a tag of type double; integer
 * fields must attach to a tag of type int.
 * Note that the use of tag data as field DOFs requires that all entities
 * with DOFs have the same number of DOFs.  That is, linear or quadratic
 * Lagrange interpolation (including p-refinement combining these two)
 * is fine.  Cubic Lagrange is not, because a vertex and an edge will
 * both have DOFs but different numbers of them, which the tag mechanism
 * can't accommodate in a single tag handle.  Usage contrary to this
 * will produce an error without affecting current field dofs or
 * attaching the field to tag data.
 * Restrictions:
 *  This function can be used for simple tensor fields or compound
 *  fields, but not for individual tensor fields that are part of a
 *  compound field.
 ******************************************************************************/
void iField_attachToTagData(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle field,
        /**< [in] Field to attach to tag data. */
    iMesh_Instance mesh_instance,
        /**< [in] iMesh instance owning the tag; must be consistent with the
        field's domain. */
    iBase_TagHandle tag_handle,
        /**< [in] iMesh tag handle containing data. */
    double* untagged_values,
        /**< [in] Value(s) to use for entities not yet tagged. This is a
        pointer because not all data is scalar. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOFStorage
 * \brief Indicate that the given field should use an external,
 *  application-created array to store DOF data.
 *
 * This function explicitly assumes that the external array data will be
 * used as the field DOFs; any existing DOFs will be destroyed by this
 * call in favor of retaining the values stored in the array. 
 * A (simple or compound) field that is stored as an array can not respond
 * to suggestions for per-entity storage; changes between block,
 * interleaved, and mixed are permitted.
 * Note that the use of array data as field DOFs requires that all
 * entities with DOFs have the same number of DOFs.  That is, linear or
 * quadratic Lagrange interpolation (including p-refinement combining
 * these two) is fine.  Cubic Lagrange is not, because a vertex and an
 * edge will both have DOFs but different numbers of them, which the tag
 * mechanism can't accommodate in a single tag handle.  Usage contrary
 * to this will produce an error without affecting current field dofs or
 * attaching the field to array data.
 * Restrictions:
 *  This function can be used for simple tensor fields or compound
 *  fields, but not for individual tensor fields that are part of a
 *  compound field.
 ******************************************************************************/
void iField_attachToArrayData(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle field,
        /**< [in] Field to attach to array data. */
    double* dof_array,
        /**< [in] Pointer to data array */
    int storage_order,
        /**< [in] Initial storage layout of the array (blocked, interleaved or
        mixed). */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOFStorage
 * \brief Give the implementation a hint about how to efficiently
 *  store data for this field.
 *
 * An implementation may or may not choose to take advantage of the hint
 * give here.  Either way, the implementation is still responsible for
 * correctly responding to all iField queries.
 * See above for description of the possible values of hint.
 ******************************************************************************/
void iField_suggestStorage(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle field,
        /**< [in] Field for which app is suggesting storage layout (ssl). */
    iField_StorageHint hint,
        /**< [in] An enumerated value for the ssl. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DFunc
 * \ingroup iField_Evaluate
 * \brief Return the values of all dfuncs at a given point.
 *
 * For a tensor_handle (and therefore distribution function kernel),
 * find all distribution function values at a particular 
 * location.  This function provides a way to get the results of
 * multiple calls to an iField_DFuncKernel in one API call.
 * dfunc_values* is the usual ITAPS array syntax.
 ******************************************************************************/
void iField_evaluateAllDF(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Tensor handle to act on. */
    const iBase_EntityHandle entity_handle,
        /**< [in] Entity in which dfuncs will be evaluated.  */
    const double coords[3],
        /**< [in] Coordinates at which to perform the evaluation. */
    const int is_param,
        /**< [in] True if coords are dfunc parametric coords; false for tensor
        physical coord system. */
    double** dfunc_values,
        /**< [in,out] Data is returned here. */
    int* dfunc_values_allocated,
        /**< [in,out] Allocated size of data array. */
    int* dfunc_values_size,
        /**< [out] Used size of data array. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DFunc
 * \ingroup iField_Evaluate
 * \brief Return the gradients of all dfuncs at a given point.
 *
 * For a tensor_handle (and therefore distribution function kernel),
 * find all distribution function gradients at a particular parametric
 * location.  This function provides a way to get the results of
 * multiple calls to an iField_DFuncKernel in one API call.
 * dfunc_values* is the usual ITAPS array syntax.
 ******************************************************************************/
void iField_evaluateGradAllDF(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Tensor handle to act on. */
    const iBase_EntityHandle entity_handle,
        /**< [in] Entity in which dfuncs will be evaluated.  */
    const double coords[3],
        /**< [in] Coordinates at which to perform the evaluation. */
    const int is_param,
        /**< [in] True if coords are parametric; false for tensor physical
        coord system. */
    const iBase_StorageOrder storage_order,
        /**< [in] Should data be blocked or interleaved? */
    double** dfunc_values,
        /**< [in,out] Data is returned here. */
    int* dfunc_values_allocated,
        /**< [in,out] Allocated size of data array. */
    int* dfunc_values_size,
        /**< [out] Used size of data array. */
    double** dfunc_grad_values,
        /**< [in,out] Data is returned here. */
    int* dfunc_grad_values_allocated,
        /**< [in,out] Allocated size of data array. */
    int* dfunc_grad_values_size,
        /**< [out] Used size of data array. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DFunc
 * \ingroup iField_Evaluate
 * \brief Return the hessians of all dfuncs at a given point.
 *
 * For a tensor_handle (and therefore distribution function kernel),
 * find all distribution function hessians at a particular parametric
 * location.  This function provides a way to get the results of
 * multiple calls to an iField_DFuncKernel in one API call.
 * dfunc_values* is the usual ITAPS array syntax.
 ******************************************************************************/
void iField_evaluateHessianAllDF(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Tensor handle to act on. */
    const iBase_EntityHandle entity_handle,
        /**< [in] Entity in which dfuncs will be evaluated.  */
    const double coords[3],
        /**< [in] Coordinates at which to perform the evaluation. */
    const int is_param,
        /**< [in] True if coords are parametric; false for tensor physical
        coord system. */
    const iBase_StorageOrder storage_order,
        /**< [in] Should data be blocked or interleaved? */
    double** dfunc_values,
        /**< [in,out] Data is returned here. */
    int* dfunc_values_allocated,
        /**< [in,out] Allocated size of data array. */
    int* dfunc_values_size,
        /**< [out] Used size of data array. */
    double** dfunc_grad_values,
        /**< [in,out] Data is returned here. */
    int* dfunc_grad_values_allocated,
        /**< [in,out] Allocated size of data array. */
    int* dfunc_grad_values_size,
        /**< [out] Used size of data array. */
    double** dfunc_hess_values,
        /**< [in,out] Data is returned here. */
    int* dfunc_hess_values_allocated,
        /**< [in,out] Allocated size of data array. */
    int* dfunc_hess_values_size,
        /**< [out] Used size of data array. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_FieldParameters
 * \brief What type of coordinate system is used for physical coordinates?
 *
 ******************************************************************************/
void iField_getCoordType(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Tensor handle to act on. */
    int* coord_type,
        /**< [out] Type of coord system. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_FieldParameters
 * \brief Retrieve name for this tensor.
 *
 ******************************************************************************/
void iField_getTensorName(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Tensor handle to act on. */
    char* name,
        /**< [in] Name of tensor; must be preallocated. */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    int* name_len
        /**< [out] Length of name string. */
);

/***************************************************************************//**
 * \ingroup iField_FieldParameters
 * \brief Retrieve distribution function kernel for this tensor.
 *
 ******************************************************************************/
void iField_getDistFunc(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Tensor handle to act on. */
    iField_DFuncKernel* dfunc,
        /**< [in] Distribution function handle for this tensor. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_TensorTemplates
 * \ingroup iField_FieldParameters
 * \brief Retrieve tensor template for this tensor.
 *
 * Restriction: tensor_handle must not be a compound field.
 ******************************************************************************/
void iField_getTensorTemplate(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Tensor handle to act on. */
    iField_TensorTemplateHandle* template_handle,
        /**< [in] Tensor template handle for this tensor. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_FieldParameters
 * \brief Retrieve units for this tensor.
 *
 * Restriction: tensor_handle must not be a compound field.
 ******************************************************************************/
void iField_getUnits(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Tensor handle to act on. */
    iField_UnitsHandle* units,
        /**< [in] Units handle for this tensor. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_FieldParameters
 * \brief Retrieve precision for this tensor.
 *
 ******************************************************************************/
void iField_getPrecision(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Tensor handle to act on. */
    int* precision,
        /**< [in] Precision for this tensor. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOF
 * \brief Returns global indices for all DOFs that affect the TensorField
 *  on this entity.
 *
 * The entity given must be one over which the dfuncs are specified.
 * For TensorFields with the same DOF associations, the order of DOFs
 * for different TensorFields should be the same (principle of least
 * astonishment).
 ******************************************************************************/
void iField_getDOFIndicesByEnt(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] Field to act on. */
    iBase_EntityHandle entity_handle,
        /**< [in] Entity for which DOF indices are desired. */
    int** indices,
        /**< [in,out] Data is returned here. */
    int* indices_allocated,
        /**< [in,out] Allocated size of data array. */
    int* indices_size,
        /**< [out] Used size of data array. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOF
 * \brief Returns global indices for all DOFs that affect the TensorField
 *  on these entities.
 *
 * The entities given must be ones over which the dfuncs are specified.
 * For TensorFields with the same DOF associations, the order of DOFs
 * for different TensorFields should be the same (principle of least
 * astonishment).
 ******************************************************************************/
void iField_getDOFIndicesByEntArr(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] Tensor to act on. */
    iBase_EntityHandle* entity_handles,
        /**< [in] Entities for which DOF indices are desired. */
    int num_entities,
        /**< [in] Number of entities passed in. */
    int** indices,
        /**< [in,out] Data is returned here. */
    int* indices_allocated,
        /**< [in,out] Allocated size of data array. */
    int* indices_size,
        /**< [out] Used size of data array. */
    int** offsets,
        /**< [in,out] Offsets for data by entity are returned here */
    int* offsets_allocated,
        /**< [in,out] Allocated size of offset array. */
    int* offsets_size,
        /**< [out] Used size of offset array. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_TensorMembers
 * \brief Destroy a tensor field.
 *
 ******************************************************************************/
void iField_destroyTensorField(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle handle,
        /**< [in] Tensor field to destroy. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Create a domain from an entity set and return its handle.
 *
 * The entities used as support for distribution functions are mesh
 * entities in the given mesh instance, and the entity set is a set of
 * mesh entities from that instance. The set may be the root set.
 ******************************************************************************/
void iField_createMeshDomainFromMeshSet(
    const iField_Instance field_instance,
        /**< [in] Field interface instance */
    const iMesh_Instance mesh_instance,
        /**< [in] iMesh instance to which the set belongs. */
    const iBase_EntitySetHandle set_handle,
        /**< [in] iMesh entity set containing the entities over which the field
        is defined. */
    iField_DomainHandle* domain_handle,
        /**< [out] Handle for the new domain. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Create a domain from a geometric entity and return its handle.
 *
 * The entities used as support for distribution functions are mesh
 * entities in the iMesh instance related to the iGeom instance by the
 * given iRel_PairHandle.  The mesh entities are related to the given
 * geometric entity, either directly or through a mesh set related to
 * the geometry entity.
 * The call sequence here is cumbersome, because both the geometric entity
 * and the iRel relation must have instances to go with them.  As far as I
 * can tell, iRel can't provide this info currently.
 ******************************************************************************/
void iField_createMeshDomainFromGeomEnt(
    const iField_Instance field_instance,
        /**< [in] Field interface instance */
    const iGeom_Instance igeom_instance,
        /**< [in] iGeom instance to which the geometric entity belongs.  */
    const iBase_EntityHandle geom_ent_handle,
        /**< [in] iGeom entity to which mesh entities are related.  */
    const iRel_Instance irel_instance,
        /**< [in] iRel instance to which the relation belongs. */
    const iRel_PairHandle relation_handle,
        /**< [in] iRel relation to the mesh entities over which the field is
        defined. */
    iField_DomainHandle* domain_handle,
        /**< [out] Handle for the new domain. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Create a domain from a set of geometric entities and return its handle.
 *
 * The entities used as support for distribution functions are mesh
 * entities in the iMesh instance related to the iGeom instance by the
 * given iRel_PairHandle.  The mesh entities are related to the geometry
 * entities in the given set of geometry entities, either directly or
 * through mesh sets related to the geometry entities.
 * The call sequence here is cumbersome, because both the geometric
 * entity set and the iRel relation must have instances to go with them.
 * As far as I can tell, iRel can't provide this info currently.
 ******************************************************************************/
void iField_createMeshDomainFromGeomSet(
    const iField_Instance field_instance,
        /**< [in] Field interface instance */
    const iGeom_Instance igeom_instance,
        /**< [in] iGeom instance to which the geometric entity belongs.  */
    const iBase_EntitySetHandle geom_set_handle,
        /**< [in] iGeom entity set to which mesh entities are related.  */
    const iRel_Instance irel_instance,
        /**< [in] iRel instance to which the relation belongs. */
    const iRel_PairHandle relation_handle,
        /**< [in] iRel relation to the mesh entities over which the field is
        defined. */
    iField_DomainHandle* domain_handle,
        /**< [out] Handle for the new domain. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Create a domain from a set of geometric entities and return its handle.
 *
 * The entities used as support for distribution functions are iGeom
 * entities in the given iGeom instance, and the entity set is a set of
 * iGeom entities from that instance. The set may be the root set.
 ******************************************************************************/
void iField_createGeomDomainFromGeomSet(
    const iField_Instance field_instance,
        /**< [in] Field interface instance */
    const iGeom_Instance geom_instance,
        /**< [in] iGeom instance to which the set belongs. */
    const iBase_EntitySetHandle set_handle,
        /**< [in] iGeom entity set containing the entities over which the field
        is defined. */
    iField_DomainHandle* domain_handle,
        /**< [out] Handle for the new domain. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Which mesh instance contains entities supporting this field?
 *
 * This function will return an error if the domain was created using
 * createGeomDomainFromGeomSet.
 ******************************************************************************/
void iField_getDomainMesh(
    const iField_Instance field_instance,
        /**< [in] Field interface instance. */
    const iField_DomainHandle domain_handle,
        /**< [in] Domain being queried. */
    iMesh_Instance* mesh_instance,
        /**< [out] Mesh instance to which mesh entities belong.  */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Which mesh set contains entities supporting this field?
 *
 * This function will return an error if the domain was created using
 * createGeomDomainFromGeomSet or createMeshDomainFromGeomSet.  If
 * createMeshDomainFromGeomEnt was used, then an error is returned if
 * the iGeom entity was associated directly to multiple mesh entities
 * rather than to an iMesh set of entities.
 * In any event, if the function returns iBase_SUCCESS, the set handle
 * must be the one which contains exactly the entities in the domain.
 ******************************************************************************/
void iField_getDomainMeshSet(
    const iField_Instance field_instance,
        /**< [in] Field interface instance. */
    const iField_DomainHandle domain_handle,
        /**< [in] Domain being queried. */
    iBase_EntitySetHandle* mesh_set_handle,
        /**< [out] Mesh set to which mesh entities belong. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Which geom ent is related to mesh entities supporting this field?
 *
 * This function will return an error if the domain was created using
 * createGeomDomainFromGeomSet, createMeshDomainFromMeshSet, or 
 * createMeshDomainFromGeomSet.
 ******************************************************************************/
void iField_getDomainGeomEnt(
    const iField_Instance field_instance,
        /**< [in] Field interface instance. */
    const iField_DomainHandle domain_handle,
        /**< [in] Domain being queried. */
    iGeom_Instance* geom_instance,
        /**< [out] iGeom instance to which geometric entity belongs.  */
    iBase_EntityHandle* geom_entity,
        /**< [out] iGeom entity related to mesh entities */
    iRel_Instance* rel_instance,
        /**< [out] iRel instance to which rel_pair belongs. */
    iRel_PairHandle* rel_pair,
        /**< [out] iRel pair handle relating geom to mesh */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Which geom set is related to mesh entities supporting this field?
 *
 * This function will return an error if the domain was created using
 * createGeomDomainFromGeomSet, createMeshDomainFromMeshSet, or 
 * createMeshDomainFromGeomEnt.
 ******************************************************************************/
void iField_getDomainGeomSet(
    const iField_Instance field_instance,
        /**< [in] Field interface instance. */
    const iField_DomainHandle domain_handle,
        /**< [in] Domain being queried. */
    iGeom_Instance* geom_instance,
        /**< [out] iGeom instance to which geometric entity belongs.  */
    iBase_EntitySetHandle* geom_entity_set,
        /**< [out] iGeom entity set related to mesh entities */
    iRel_Instance* rel_instance,
        /**< [out] iRel instance to which rel_pair belongs. */
    iRel_PairHandle* rel_pair,
        /**< [out] iRel pair handle relating geom to mesh */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Which geom set contains geom entities supporting this field?
 *
 * This function will return an error if the domain was created using
 * createMeshDomainFromGeomSet, createMeshDomainFromMeshSet, or 
 * createMeshDomainFromGeomEnt.
 ******************************************************************************/
void iField_getGeomDomainInfo(
    const iField_Instance field_instance,
        /**< [in] Field interface instance. */
    const iField_DomainHandle domain_handle,
        /**< [in] Domain being queried. */
    iGeom_Instance* geom_instance,
        /**< [out] iGeom instance to which geometric entity belongs.  */
    iBase_EntitySetHandle* geom_entity_set,
        /**< [out] iGeom entity set related to mesh entities */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_EntityIterator
 * \brief Initialize an iterator over specified entity type and topology
 *
 * Initialize an array iterator over specified entity type and topology.
 * Iterator returned can be used as input to functions returning entities 
 * for the iterator.  If all entities of a specified type and/or topology
 * are to be iterated, specify iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, 
 * respectively.  Specified type or topology must be a value in the
 * iBase_EntityType or iMesh_EntityTopology enumerations, respectively.
 * Note that we will probably want (certainly need?) different syntax
 * for iterators over geometric domains.
 ******************************************************************************/
void iField_initEntIter(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_DomainHandle domain_handle,
        /**< [in] Domain to iterate over */
    const int requested_entity_type,
        /**< [in] Type of entity to iterate */
    const int requested_entity_topology,
        /**< [in] Topology of entity to iterate */
    iMesh_EntityIterator* entity_iterator,
        /**< [out] New iterator returned by function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_EntityIterator
 * \brief Initialize an array iterator over specified entity type, topology, and 
 *  size
 *
 * Initialize an array iterator over specified entity type, topology, and 
 * size, for a specified domain.  Iterator returned can be used 
 * as input to functions returning entities for the iterator.  If all 
 * entities of a specified type and/or topology are to be iterated, 
 * specify iBase_ALL_TYPES or iMesh_ALL_TOPOLOGIES, respectively.  
 * Specified type or topology must be a value in the iBase_EntityType or 
 * iMesh_EntityTopology enumerations, respectively.
 * Note that we will probably want (certainly need?) different syntax
 * for iterators over geometric domains.
 ******************************************************************************/
void iField_initEntArrIter(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_DomainHandle domain_handle,
        /**< [in] Domain to iterate over */
    const int requested_entity_type,
        /**< [in] Type of entity to iterate */
    const int requested_entity_topology,
        /**< [in] Topology of entity to iterate */
    const int requested_array_size,
        /**< [in] Size of blocks for iterator to return.  */
    iMesh_EntityArrIterator* entArr_iterator,
        /**< [out] New iterator returned by function */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Destroy a domain.
 *
 * This should refuse to destroy if the domain is in use.  This only
 * destroys the Domain, without affecting the underlying mesh
 * entities, or any possible mesh sets, geom ents, or geom sets.
 ******************************************************************************/
void iField_destroyDomain(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_DomainHandle handle,
        /**< [in] Domain to destroy. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_CoordinateField
 * \brief Create a coordinate field and return its handle.
 *
 * Coordinate fields are different from other fields:  we know in
 * advance that they are vector fields measuring length.  This
 * simplifies the calling sequence significantly
 * Note that coordinate fields can -only- be evaluated in parametric
 * coordinates.  For one thing, otherwise you're asking "What are the
 * values of (x,y,z) for this value of (x,y,z)?"  Pretty sure I know the
 * answer to that question....
 * MCM: We don't want to limit ourselves to length. We ought to be
 * able to define fields on any space including fields on pressure
 * temperature, 2D phase space for a material.
 * CFOG reply:  I don't disagree, but let's make the simple things
 * simple and add that as a separate function later.
 ******************************************************************************/
void iField_createCoordinateField(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iMesh_Instance mesh_instance,
        /**< [in] The mesh this coordinate field will provide coordinates for.  */
    const iField_UnitsHandle units,
        /**< [in] Units to measure length in. */
    int coord_type,
        /**< [in] Coordinate system to use. */
    const iField_DFuncKernel dfunc_handle,
        /**< [in] Distribution function kernel to use for coordinates.  */
    iField_Handle* field_handle,
        /**< [out] Newly create coordinate field. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_TensorMembers
 * \brief Return the component tensor fields making up a compound field.
 *
 ******************************************************************************/
void iField_getTensorFields(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    const iField_TensorHandle tensor_handle,
        /**< [in] Field to operate on. */
    iField_TensorHandle** tensor_handles,
        /**< [in,out] Array of component tensor fields. */
    int* tensor_handles_allocated,
        /**< [in,out] Allocated size of TF array. */
    int* tensor_handles_size,
        /**< [out] Used size of TF array. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOF
 * \brief Retrieve DOFs associated with this entity.
 *
 * This function returns dofs associated directly with this entity and,
 * optionally, entities on its closure.  This distinction is needed to
 * support mixed-accuracy discretizations, although full support for
 * this is not yet in the interface.
 * For collective calls, DOFs are returned in either blocked,
 * interleaved, or mixed order, as requested.
 * NOTE: With this function signature, the implementation is responsible
 * for returning a pointer to the location of that data; the application
 * is responsible for correctly interpreting the data at that address.
 * If (a pointer to) iField_DEFAULT_PRECISION is passed in, the
 * implementation will return data at the highest precision used by any
 * component field tensor and set *precision to this value.
 * Array sizes are given in terms of number of values of the field's
 * data type, at *precision, NOT in terms of bytes.  This implies that
 * pointer conversion will make pointer arithmetic easier / more
 * readable. 
 * NOTE: The syntax below is ridiculously awkward to work with in
 * practice.  In the near future, there will be precision-specific
 * versions of all get/setDOF and field evaluation functions, so that
 * conversion between data pointed to by a void* and double (or
 * whatever) actual data.  Implementation of this should be possible
 * relatively easy using pre-processor tricks to do the equivalent of
 * C++ templates.
 * This approach will not only simplify app programming, but also allow
 * compile time type checking (no passing float data incorrectly labeled
 * as double....).
 * We will also REQUIRE that if a type-specific evaluation method is
 * implemented then it is implemented in the precision the interface
 * implies.
 * MCM: We could defer the issue of low-level datatype variability and
 * put in as a place-holder a simple macro DTYPE() which expands to
 * two things, an enum indicating type and a void* pointer but treat
 * that all as double for initial implementations.
 ******************************************************************************/
void iField_getDOFbyEnt(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] The tensor being read. */
    const iBase_EntityHandle entity_handle,
        /**< [in] Entity for which DOFs are requested. */
    const int include_closure,
        /**< [in] DOFs for closure, too? */
    int* precision,
        /**< [out] Precision of returned data. */
    void** dofs,
        /**< [in,out] Storage for dof data. */
    int* dofs_allocated,
        /**< [in,out] Allocated size, of values of size precision. */
    int* dofs_size,
        /**< [out] Used size, of values of size precision. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOF
 * \brief Retrieve DOFs associated with an array of entities
 *
 ******************************************************************************/
void iField_getDOFbyEntArr(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] The tensor being read */
    const iBase_EntityHandle* entity_handles,
        /**< [in] Entities for which DOFs are requested */
    const int num_handles,
        /**< [in] number of entities in entity_handles array */
    const int include_closure,
        /**< [in] DOFs for closure too? */
    const int storage_order,
        /**< [in] desired storage order of returned data */
    int* precision,
        /**< [out] precision of returned data */
    void** dofs,
        /**< [in,out] Array of returned DOF data */
    int* dofs_allocated,
        /**< [in,out] allocated size of dofs array */
    int* dofs_size,
        /**< [out] occupied size of dofs array */
    int** offsets,
        /**< [in,out] array of offsets into DOFs array for each entity's DOFs */
    int* offsets_allocated,
        /**< [in,out] allocated size of offsets array */
    int* offsets_size,
        /**< [out] occupied size of offsets array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOF
 * \brief Retrieve DOFs associated with all entities in the domain 
 *
 ******************************************************************************/
void iField_getAllDOF(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] The tensor being read */
    int* precision,
        /**< [out] precision of returned data */
    const int storage_order,
        /**< [in] desired storage order of returned data */
    void** dofs,
        /**< [in,out] Array of returned DOF data */
    int* dofs_allocated,
        /**< [in,out] allocated size of dofs array */
    int* dofs_size,
        /**< [out] occupied size of dofs array */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOF
 * \brief Set DOFs associated with this entity.
 *
 * This function stores dofs associated directly with this entity and,
 * optionally, entities on its closure.  This distinction is needed to
 * support mixed-accuracy discretizations, although full support for
 * this is not yet in the interface.
 * For collective calls, DOFs are given in either blocked, interleaved,
 * or mixed order, as requested.
 * NOTE: The application is responsible for sending data at the
 * specified precision; the implementation is responsible for correctly
 * interpreting the data and storing it at the precision specified at
 * field creation.
 * Array sizes are given in terms of number of values at *precision, NOT
 * in terms of bytes.  
 ******************************************************************************/
void iField_setDOFbyEnt(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] The tensor being written. */
    const iBase_EntityHandle entity_handle,
        /**< [in] Entity for which DOFs are requested. */
    const int include_closure,
        /**< [in] DOFs for closure, too? */
    const int storage_order,
        /**< [in] storage order of dof array */
    int precision,
        /**< [in] Precision of input data. */
    void* dofs,
        /**< [in] Storage for dof data. */
    const int dofs_size,
        /**< [in] Used size, in bytes. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOF
 * \brief Brief unavailable
 *
 * Description unavailable
 ******************************************************************************/
void iField_setDOFbyEntArr(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< description unknown */
    const iBase_EntityHandle* entity_handles,
        /**< description unknown */
    const int num_entities,
        /**< description unknown */
    const int include_closure,
        /**< description unknown */
    const int storage_order,
        /**< description unknown */
    int precision,
        /**< description unknown */
    void* dofs,
        /**< description unknown */
    const int dofs_size,
        /**< description unknown */
    int* offsets,
        /**< description unknown */
    const int offsets_size,
        /**< description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DOF
 * \brief Brief unavailable
 *
 * Description unavailable
 ******************************************************************************/
void iField_setAllDOF(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< description unknown */
    int precision,
        /**< description unknown */
    const int storage_order,
        /**< description unknown */
    void* dofs,
        /**< description unknown */
    int* dofs_size,
        /**< description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Evaluate
 * \brief Evaluate a compound field at a point in a known entity.
 *
 * Knowing which entity should allow significant optimization.
 * In the collective call, each set of coordinates must have an entity
 * handle to go with it.
 * NOTE: This function returns all data at the same precision.  The
 * implementation is responsible for returning a pointer to the location
 * of that data; the application is responsible for correctly
 * interpreting the data at that address.  If (a pointer to)
 * iField_DEFAULT_PRECISION is passed in, the implementation will return
 * data at the highest precision used by any component field tensor and
 * set *precision to this value.
 * Array sizes are given in terms of number of floating point numbers at
 * *precision, NOT in terms of bytes.  This implies conversion to a
 * pointer to floating point number before doing pointer arithmetic.
 * CFOG Oct 4, 2010
 *   To handle evaluation of solution, gradient, Hessian, etc for a
 *   field, do we want to add functions, or add an arg to the current
 *   functions specifying what to evaluate?
 ******************************************************************************/
void iField_evaluateInEnt(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] Tensor instance to act on. */
    const iBase_EntityHandle entity_handle,
        /**< [in] Entity in which dfuncs will be evaluated.  */
    const double* coords,
        /**< [in] Coordinates at which to evaluate dfunc.  */
    const int is_param,
        /**< [in] True if coords are parametric; false for tensor physical
        coord system. */
    int* precision,
        /**< [in,out] Precision of output data. */
    void** field_values,
        /**< [in,out] Data is returned here. */
    int* field_values_allocated,
        /**< [in,out] Allocated size of data array. */
    int* field_values_size,
        /**< [out] Used size of data array. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Evaluate
 * \brief Brief unavailable
 *
 * Description unavailable
 ******************************************************************************/
void iField_evaluateInEntArr(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< description unknown */
    const int num_points,
        /**< description unknown */
    const iBase_EntityHandle* entity_handles,
        /**< description unknown */
    const double* coords,
        /**< description unknown */
    const int is_param,
        /**< description unknown */
    int* precision,
        /**< description unknown */
    void** field_values,
        /**< description unknown */
    int* field_values_allocated,
        /**< description unknown */
    int* field_values_size,
        /**< description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Evaluate
 * \brief Evaluate a compound field at a point near a known entity.
 *
 * This is the case that will be intermediate in efficiency.
 * Note that, for this case, coordinates are in the physical coordinate
 * frame: if you're using parametric coordinates for a known entity, use
 * evaluateInEnt even if the point lies physically outside the entity.
 * See notes for evaluateInEnt
 ******************************************************************************/
void iField_evaluateNearEnt(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] Tensor to act on. */
    const iBase_EntityHandle entity_handle,
        /**< [in] Entity near which dfuncs will be evaluated.  */
    const double* coords,
        /**< [in] Coordinates at which to evaluate dfunc.  */
    int* precision,
        /**< [in,out] Precision of output data. */
    void** field_values,
        /**< [in,out] Data is returned here. */
    int* field_values_allocated,
        /**< [in,out] Allocated size of data array. */
    int* field_values_size,
        /**< [out] Used size of data array. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Evaluate
 * \brief Brief unavailable
 *
 * Description unavailable
 ******************************************************************************/
void iField_evaluateNearEntArr(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< description unknown */
    const int num_points,
        /**< description unknown */
    const iBase_EntityHandle* entity_handles,
        /**< description unknown */
    const double* coords,
        /**< description unknown */
    const int is_param,
        /**< description unknown */
    int* precision,
        /**< description unknown */
    void** field_values,
        /**< description unknown */
    int* field_values_allocated,
        /**< description unknown */
    int* field_values_size,
        /**< description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Evaluate
 * \brief Evaluate a compound field with no location hint.
 *
 * This is the expensive case, because we have to find the entity first.
 * Might as well return it for future use....
 * Note that, for this case, coordinates are in the physical coordinate
 * frame: if you're using parametric coordinates for a known entity, use
 * evaluateInEnt even if the point lies physically outside the entity.
 * See notes for evaluateInEnt
 ******************************************************************************/
void iField_evaluate(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] Field to act on. */
    const double* coords,
        /**< [in] Coordinates at which to evaluate dfunc.  */
    iBase_EntityHandle* entity_handle,
        /**< [out] Entity in which dfuncs were evaluated.  */
    void** field_values,
        /**< [in,out] Data is returned here. */
    int* field_values_allocated,
        /**< [in,out] Allocated size of data array. */
    int* field_values_size,
        /**< [out] Used size of data array. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Evaluate
 * \brief Brief unavailable
 *
 * Description unavailable
 ******************************************************************************/
void iField_evaluateArr(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< description unknown */
    const int num_points,
        /**< description unknown */
    const double* coords,
        /**< description unknown */
    iBase_EntityHandle** entity_handles,
        /**< description unknown */
    int* entity_handles_allocated,
        /**< description unknown */
    int* entity_handles_size,
        /**< description unknown */
    void** field_values,
        /**< description unknown */
    int* field_values_allocated,
        /**< description unknown */
    int* field_values_size,
        /**< description unknown */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Domains
 * \brief Retrieve information about the domain for this compound field.
 *
 ******************************************************************************/
void iField_getDomain(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] Tensor to act on. */
    iField_DomainHandle* domain_handle,
        /**< [out] The domain for this Field. */
    int* err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_Persistence
 * \brief Save a Field to a file.
 *
 * There is a problem with maintaining the binding between iField and
 * iMesh (or iGeom) persistently. In a running executable, we've got
 * the iField and iMesh instance handles. But, when loading from
 * disk, all we've got is filenames.
 * Our options for maintaining iMesh<->iField binding persisently are...
 *    a) Store iMesh file name in iField (perhaps reverse as well).
 *       This will probably require methods specifically for this task
 *       in both interfaces.
 *    b) Store iMesh and iField instances to the same file. This is
 *       the most natural thing to do anyways. But, I don't think we
 *       can assume all implementations are prepared to handle this
 *       so I am not sure we can make this a requirement. 
 *    c) Rely upon some naming convention such that given the name
 *       of one piece, you can construct the other and vice versa.
 *       This is the simplest strategy but feels kludgy. 
 *    d) Ignore the problem. Let the applications worry about it. I
 *       think this is just asking for trouble.
 *    e) A 'meta' file containing the iField/iMesh filenames.
 *       An xml-ish meta file might serve us well here.
 * Note: the above assumes all fields are in one file. But the save
 * call defined here doesn't necessarily REQUIRE that. Perhaps it
 * should. If we somehow loose track of the mesh file, then all of
 * our field data is pretty much totally useless.
 * If this function is called multiple times for different fields in the
 * same iField_instance what happens? Can all the data go to the same
 * file? Must it go into different files? If think if the latter is true,
 * we're asking for serious I/O performance problems. On the other hand,
 * if the caller passes a different filename each time, then the impl.
 * has no choice, either.
 * CFOG: Checksums on iMesh file would be useful to keep to know if
 * something goes wrong with the mesh file. Another issue is mixing,
 * for example, a GRUMMP iMesh instance with a MOAB iField instance.
 * There is no reason this should NOT be supported. Nonetheless, this
 * probably does raise some issues regarding save/load operations that
 * we've yet to consider.
 * File format and options will inevitably be implementation dependent.
 ******************************************************************************/
void iField_save(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] tensor to act on. */
    const char* filename,
        /**< [in] What file to write to. */
    const char* options,
        /**< [in] Options, probably implementation dependent. */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    const int filename_len,
        /**< [in] Length of filename string. */
    const int options_len
        /**< [in] Length of options string. */
);

/***************************************************************************//**
 * \ingroup iField_Persistence
 * \brief Read a Field definitions from a file 
 *
 * File format and options will inevitably be implementation dependent.
 * Note: No actual field dofs are read into memory in this call.
 * MCM: maintaining the binding between dofs (in iField) and entities (in iMesh)
 * (e.g. the dof associations) involves assumptions regarding the ordering
 * of these things in the buckets in which they are stored. iField may need
 * to impose restrictions on what iMesh is allowed to do in the way of
 * re-ordering things on its _save _load operations. Or, iMesh needs to
 * deliver some permutation vector to iField so that iField can re-order its
 * dofs to match iMesh (for those dofs that are not handled as tags).
 * MCM: The above sounds like it is similar to the problem of how do we deal
 * with addition/deletion of entities in iMesh? What happens to the dofs associated
 * with those entities? I don't think we can prevail upon iMesh to inform
 * iField that the associated dofs also need to be deleted. However, we
 * might be able to require iField to ask the iMesh instance if the entities
 * for which it is storing dofs have 'changed'. But even so, knowing the
 * ordering has changed is not sufficient. iField would need to query from
 * iMesh the new ordering and then adjust dof storage to match it. I think
 * having call in which all this magic was buried would be helpful to the
 * iField client. Something like...
 * void iField_syncDofStorage(const iField_Instance instance, int *err);
 * might be sufficient.
 ******************************************************************************/
void iField_load(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle tensor_handle,
        /**< [in] Tensor to act on. */
    const char* filename,
        /**< [in] What file to read from. */
    const char* options,
        /**< [in] Options, probably implementation dependent. */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    const int filename_len,
        /**< [in] Length of filename string. */
    const int options_len
        /**< [in] Length of options string. */
);

/***************************************************************************//**
 * \ingroup iField_Persistence
 * \brief Read dofs for (some) fields in the instance
 *
 * MCM: This approach assumes an impl. can allocate memory for dofs
 * separately from the field 'header' information itself.
 * MCM: We discussed adjusting the iField_createTensorField call to
 * accept an option to indicate if the memory should be allocated
 * then or deferred to later. That is basically the same idea as
 * we are hoping to achieve here in the load call.
 ******************************************************************************/
void iField_loadDofs(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle fields_whose_dofs_to_load,
        /**< [in] The field handles for which dofs should be loaded */
    const char* options,
        /**< [in] user defined options */
    int* err,
        /**< [out] Returned Error status (see iBase_ErrorType) */
    const int options_len
        /**< [in] length of the options string */
);

/***************************************************************************//**
 * \ingroup iField_DFunc
 * \brief Specify a particular dfunc for an entity.
 *
 * This function can be called on compound fields or simple fields that
 * are not part of compound fields.
 *
 * New DOFs aren't set by this call (what values should we guess?), so
 * the app is going to have to do that after.
 ******************************************************************************/
void iField_changeDFunc(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle field,
        /**< [in] The tensor field handle being changed. */
    iBase_EntityHandle entity,
        /**< [in] The entity handle being changed. */
    iField_DFuncKernel new_dfunc,
        /**< [in] The new dfunc to set. */
    int *err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DFunc
 * \brief Specify a particular dfunc for an array of entities.
 *
 * This function can be called on compound fields or simple fields that
 * are not part of compound fields.
 *
 * New DOFs aren't set by this call (what values should we guess?), so
 * the app is going to have to do that after.
 ******************************************************************************/
void iField_changeDFuncArr(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle field,
        /**< [in] The tensor field handle being changed. */
    iBase_EntityHandle entities,
        /**< [in] The array of entities being changed. */
    const int num_entities,
        /**< [in] The size of the entities array. */
    iField_DFuncKernel new_dfunc,
        /**< [in] The new dfunc to set. */
    int *err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \ingroup iField_DFunc
 * \brief Return all dfuncs in use for any entities for this field.
 *
 ******************************************************************************/
void iField_getAllDFuncs(
    const iField_Instance instance,
        /**< [in] iField instance handle */
    iField_TensorHandle field,
        /**< [in] The tensor field handle being queiried. */
    iField_DFuncKernel **dfunc,
        /**< [in,out] The returned array of dfuncs in use. */
    int *dfunc_allocated,
        /**< [in,out] Allocated size of dfuncs array. */
    int *dfunc_size,
        /**< [out] Occupied size of dfuncs array. */
    int *err
        /**< [out] Returned Error status (see iBase_ErrorType) */
);

/***************************************************************************//**
 * \defgroup iField iField 
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_Initialization Initialization
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_ErrorHandling Error Handling
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_TensorTemplates Tensor Templates
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_TensorMembers Tensor Member Fields
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_CoordinateField Coordinate Field
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_DOFStorage DOF Storage Options
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_DFuncs Distribution Functions
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_FieldParameters Definition Paramaters of a Field
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_Domains Domain of a Field
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_EntityIterator Entity Iterators
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_DOF Setting and Getting Degrees of Freedom (DOFs)
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_Evaluate Evaluating Fields 
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**
 * \defgroup iField_Persistence Saving and Loading fields to persistent storage
 * \ingroup iField
 ******************************************************************************/

/***************************************************************************//**

\page ifield iField: ITAPS Field Interface

\section INTRO Overview

The iField interface is intended to provide support for tensor fields
and collections of tensor fields, with a strong bias towards supporting
first and best the capabilities required for scientific computing
applications.

Scientific computing applications deal in physical quantities expressed
as tensors: scalars such as temperature, vectors such as velocity, and
second-order tensors such as stress. In practice, these are formally
tensor fields: a tensor field assigns a tensor to each point in a
mathematical space (typically a Euclidean space or manifold).

Because tensors have a geometric interpretation, their underlying
physical meaning is independent of the coordinate system in which they
are defined; the numerical value of a tensor depends on the coordinate
system, and so any numerical instantiation of a tensor field is
incomplete without a specification of its coordinate system.

The ITAPS Field interface represents a field as a linear combination of
distribution functions and degrees of freedom (dof's): $\sum_{i} w_{i}
f_{i} \left(\vec{x}\right)$.  For typical discrete fields, the
distribution functions $f_{i}$ have compact support over a single
highest-dimensional mesh entity, and are each more or less associated
with a mesh entity, which may be of lower dimension; for instance, the
interpolation functions used in the finite element method may be
associated with a vertex, edge, face, or region. These distribution
functions are known in advance. The degrees of freedom $w_{i}$ are
associated with both an entity and a distribution function; calculation
of these degrees of freedom is the goal of a typical scientific computing
application.

An ITAPS Field is, for practical purposes, a function. The
representation of the function will, in most use cases, be piecewise and
will not necessarily by continuous (never mind smooth), but with
reasonable definitions in a few tight spots, a Field is in fact a
function.

In addition to supporting physical tensor fields, there's a reasonable
argument for allowing integer or boolean values for fields.  While these
don't represent physical tensors, and need to be restricted to scalars
(coordinate transforms for vectors of booleans, anyone?), applications
do and will continue to want to use such things.  This case implies
certain restrictions on distribution functions for integer and boolean
scalar "fields" and on the operations one can perform on them.

To support this representation of Fields and the sorts of operations on
Fields that applications require, we need precise definitions of the
domain of a function, of the data that a Field represents, and of how we
relate the domain and range for a Field.

\section ADM Abstract Data Model

\subsection DOM Domain

\subsubsection Preliminary Definitions

Primitive SET (pSET) An infinite point set corresponding, topologically,
  to an iMesh or iGeom Entity.

  Examples:

    iMesh_TRIANGLE
    iMesh_HEXAHEDRON

The key feature of a pSET is its topology and the adjacency information
implied by, for instance, knowing the identities of the faces, edges,
and vertices on the closure of a region.  Geometric information
(locations of vertices, shapes of higher-dimensional entities) is a
separate matter.

Pre-Domain: The union of a collection of primitive point sets; that is,
  an iMesh or iGeom EntitySet.

A pre-domain is a topological description of the domain; specifically,
it has no coordinates associated with it to define the spatial location
of iMesh and iGeom Entities. It represents an infinite point set on
which the field will be defined, decomposed into iMesh/iGeom primitives
(pSETs) with adjacency information to knit those pSETs together to cover
the infinite point set.



Domain (mesh based)

The domain of a field is expressed, at the lowest level, by a collection
of mesh entities over which the field is defined.  This collection can
potentially be expressed in multiple ways: as the entities in an iMesh
entity set (in many use cases, the root set) or as the mesh entities
related (by iRel) to one or more iGeom entities.

Loose coupling between mesh topology and geometry is a benefit in
several cases --- including mesh smoothing/relaxation, shape
optimization, and deforming domain problems.  This suggests that some
cooperation between iMesh and iField is ideal for handling
coordinates, with some cases handled strictly through iMesh, others
strictly through iField, and some jointly.

For practical reasons, the iMesh specification defines one very
important field; the coordinate field. In addition, the iMesh
specification defines only the most common coordinate field
representation; a piecewise linear, Lagrange field. Support for more
sophisticated coordinate fields will be possible but only through
iField, not through iMesh's 'native' coordinate field representation.

When a new mesh is instantiated, an iMesh client will indicate if the
mesh will employ iMesh's 'native' coordinate field representation or not
via a new boolean argument to iMesh_newMesh()

void iMesh_newMesh(
    const char *options,
        Pointer to implementation-specific options string
    int lagrange_coords,
        flag indicating if iMesh' lagrange coordinates will be used
    iMesh_Instance *instance,
        Pointer to iMesh instance handle returned from function
    int *err,
        Pointer to error type returned from function
    int options_len
        Length of the character string pointed to by options
);

If 'int lagrange_coords' is zero (false)...

    that means the coordinate field for the mesh will be handled
    through iField. In this case, the iMesh implementation is free
    to do whatever it wishes with coordinate values passed via
    createVtx calls including ignoring that data entirely and not
    even storing it. Calls to either iMesh_getVtxArrCoords or
    iMesh_getVtxCoord will return no coordinate data and set the err
    return to iBase_NO_LAGANGE_COORDS. Calls to
    iMesh_setVtxArrCoords and iMesh_setVtxCoord() will make no
    attempt to interpret coordinate data passed in and will likewise
    return err set to iBase_NO_LAGANGE_COORDS. Calls to
    iMesh_createVtxArr and iMesh_createVtx will make no attempt to
    interepret coordinate data passed in (e.g. will not attempt to
    de-reference pointers to coordinate values), but will return
    iBase_SUCCESS.
        
    Note, even if iField is later used to define a coordinate field
    that is in fact a piecewise-linear, lagrange coordinate field
    for the mesh, there will be no obligation that that fields dofs
    make its way into the iMesh mesh instance such that it be
    managed by iMesh as it otherwise ordinarily would have.

Otherwise, 'int lagrange_coords' is non-zero (true) and...

    that means the coordinate field for the mesh will be handled
    BOTH through iMesh AND through iField. By 'BOTH', here that does
    not mean that data is redundantly stored. It only means that it
    is redundantly accessible through set/getVtx methods in iMesh as
    well as set/getDof methods in iField.  For some implementations,
    this will require a different access/storage method for this
    particular type of coordinate field than for other fields.

\subsection TENSOR Representing Tensor Data

Tensor Templates

A tensor template is intended to define all of the information about a
physical tensor that is dependent only on the physical quantity being
represented.  For instance, a 3D velocity vector is always a 3-vector of
real numbers and represents a length/time.  A 2D stress tensor is always
a 2nd-order tensor (2x2) of real numbers and represents a force/area
(mass-length/(time^2 length^2)).  Electric field is a 3D vector of real
numbers and represents a force/unit charge.  These statements will be
true for any 3D velocity or 2D stress or 3D electric field, regardless
of system of units, coordinate system, precision of stored values, etc.

Important things to know about a tensor template

  Tensor order (integer) 0 for scalars, 1 for vectors, 2 for
  tensors mapping one vector to another, etc.

  Number of tensor components (integer) Is this a 2-vector, a
  3-vector, etc.

  Algebraic type (enumeration) This describes the rules for
  combining values and/or the range of legal values. For
  instance, real numbers, complex numbers, barycentric coordinate
  vectors, and quaternions all have different algebraic rules.

  Quantity (handle) See definition later in this document.

\subsection DFUNC Distribution Functions on Meshes

At this stage, the focus is on making the easy things easy. The easiest
case for distribution functions on meshes is the case where each entity
has the same set of distribution functions; that is, the case in which
the field representation is the same for every highest (topological)
dimensional entity.  Note that this does not preclude things like linear
elements for pressure and quadratic for velocity (these are, after all,
different fields and can have different dfuncs).  Nor does this preclude
things like edge-based elements: these split a typical shape function
over multiple regions (as do other continuous FE methods).  Note also
that the definition refers to the highest TOPOLOGICAL dimension entity:
this implies that a field can be defined over a collection of faces
(say, all or part of a boundary), with dfuncs defined over the faces.

Distribution functions and their derivatives will be hard-coded in a
distribution function kernel (DFK), which will contain all information
required to evaluate all distribution functions (and their gradients,
etc) for that highest-dimensional entity (HDE).  In practice, this takes
the form of a callback from the implementation to some external
distribution function kernel, which will evaluate a specified dfunc (one
of perhaps many for the HDE) as a specified location in a specified
entity.

In this scenario, each distribution function (dfunc) has support over
exactly one highest (topological) dimensional entity; one or more
distribution functions (dfuncs) can be have the same highest-dimensional
entity as their support. Each distribution function has exactly one
degree of freedom (dof) associated with it per tensor component. The
dfunc and dof may be associated directly with a vertex, edge, or face
instead of a region (think high-order FE). The dofs and dfuncs must have
a common ordering, dictated by the DFK and understood by the
application.  The implementation needs only to know how to call the DFK
callback and sum dof*dfunc value at a point.

As a practical matter, it is obviously very important to know how to
convert global coordinates into local coordinates and vice versa, as the
dfuncs will necessarily have to use local coordinates. Obvious choices
like barycentric and shifted cartesian coordinates could be specified in
an enum, but the current API anticipates that there will be a wide and
weird enough range of local coordinate systems that it makes sense for
the DFK to include functions for coordinate conversion in both
directions, as well as the Jacobian of the transformation.

The More General Case

More generally, it's often desirable to define distribution functions
that are not associated with highest-dimensional entities, or that are
have a support that extends beyond a single HDE.  Examples of this include
vertex-centered finite-volume schemes (dfunc associated with a vertex,
with support over parts of all regions incident on that vertex),
spectral-like methods, and filtering or convolution kernels.  Although
details have not yet been worked out, the interface for supporting these
types of DFK's is likely to bear some syntactic and semantic resemblence
to the iMeshP ghost descriptions.

Examples (all dof counts for tets at this point...)

  Galerkin FE (quadratic)
     One dof/vertex, one dof/edge, zero dof/face, zero dof/region.
     Local coordinates: barycentric.
     Code for 10 dfuncs (standard shape functions) must be provided.

  Discontinuous Galerkin (cubic)
     Zero dof/vertex, zero dof/edge, zero dof/face, twenty dof/region.
     Local coordinates: barycentric.
     Code for 20 dfuncs (standard shape functions) must be provided.

  Cell-centered FV (quadratic)
     Zero dof/vertex, zero dof/edge, zero dof/face, 10 dof/region.
     Local coordinates: cartesian w/ origin at centroid.
     Code for 10 dfuncs (monomials) must be provided.

  Vertex-centered FV (quadratic)
     10 dof/vertex, zero dof/edge, zero dof/face, zero dof/region.
     Local coordinates: cartesian w/ origin at vertex(?), plus a filter
       function that selects which vertex-centered control volume that
       covers part of a region the point actually falls into (filter
       function uses local barycentric coords).
     Code for 10 dfuncs (monomials) must be provided.
     Note: this is a bit of a hack; in practice most field evals that a
       solver does are already pre-identified with a CV.  This is a case
       that cries out for a more sophisticated approach to defining
       support of dfuncs.

     
\subsection Tensor Field

A tensor field is an intermediate step in instantiating a tensor
template: at this level, units, data precision, a type for physical
coordinates, a distribution function kernel and a domain are added.

  Units (handle) See definition later in this document.

  Coordinate type (enum): Possible values include cartesian,
  cylindrical, and spherical. This is critical for tensor
  transformations.

  DFK   See above.

  Data precision.  This identifies whether floating point data is stored
  in single, double, or quad precision.  This enumeration may need to be
  expanded.

  Domain


In addition to fields that describe single physical tensors, it is often
convenient in scientific computing to have multiple tensor fields stored
and accessed together.  For instance, in incompressible flow, a vector
velocity and scalar pressure are both required, and often stored
interleaved (uvwPuvwP...) to improve cache hit rates.  To support this
common usage --- including making it relatively easy for implementations
to store variables together --- the iField API includes a function to
create a collection of tensor fields together as well as an overall
field that acts as a container for these multiple tensors.  In a
possible abuse of the data model, these compound fields use the same
data type as regular tensor fields; implementations are required to be
able to distinguish between the two.  Some queries for a tensor field
can not be reasonably answered for a compound field; these cases are
noted in the API spec.

Glossary

Coordinate system How space is laid out and measured. A
  coordinate system may be local or global, and either may be
  physical or parametric. Common examples will surely include
  cartesian, polar, spherical, and barycentric coordinates.
  User-defined coordinate systems can be infinite in variety,
  obviously. Each coordinate field will be associated with a
  coordinate system; a coordinate system can apply to many
  coordinate fields. More on this when coordinate fields come
  into the interface in earnest.

  Key functionality: conversion to/from cartesian, plus producing
  the metric and Jacobian info required for differentiation and
  integration.

Quantity Measurable physical property of phenomena/bodies/substances.
	 
    Examples: mass, length, time, acceleration, volume, power,
      force, illuminance

  Any Quantity can be expressed in terms of the seven basic physical
  quantities, corresponding to the seven base units of the SI system:
  mass, length, time, temperature, electric charge, luminous intensity
  and amount of substance. A Quantity has the following data

  (required) int[7] num_powers; Powers of 7 basic quantities in
    numerator

  (required) int[7] den_powers; Powers of 7 basic quantities in
    denominator

  (required) int flags; bool bits to indicate e.g.
    intensive/extensive, conserved, etc.

  (optional) char *name; name of quantity

  (optional) char *abbr; abbreviation of quantity.

  Operations: 7 basic quantities are pre-defined and guaranteed
  to always exist. New, derived, quantities are defined by
  multiplying together old ones.

  Notes:

  It useful to distinguish a quantity representing a mass
    fraction (e.g. mass/mass) from a quantity representing a
    volume fraction (e.g. volume/volume) for example.

  Essentially there are two kinds of quantities; base quantities
  and derived quantities. Base quantities are, by necessity,
  defined by the data model specification itself and are
  guaranteed to exist. Derived quantities shall be definable by
  applications as necessary. See
  physics.nist.gov/cuu/Units/units.html for detailed discussion.

  Unit An adopted convention for assigning a numerical value to a
  Quantity. A Unit has

  (required) handle quant; the Quantity the unit is a measure for

  (required) double scale; scale of this unit relative to the
    default unit of measure for quant.

  (optional) double offset; likewise for offset units

  (optional) double logbase; likewise for log-related units

  (optional) double logcoef; likewise for log-related units

  (optional) char *name;

  (optional) char *abbr;

  Examples:

  slug is a unit of measure of the quantity mass

  pound is a unit of measure of the quantity force

  minute is a unit of measure of the quantity time

  meters/second is a unit of measure of the quantity velocity

  furlongs/fortnight is a unit of measure of the quantity
    velocity

  decibels is a (log scaled) unit of measure of the quantity
    power.

  degrees Fahrenheit is a (offset) unit of measure of the
    quantity thermodynamic temperature.

  Operations:

  Associate a Unit (of measure) with a Quantity

  Combine Unit A and Unit B producing new Unit C by scale factor
    and optional power.

  Set log and/or offset scaling of Unit A relative to Unit B

  Unit System A collection of quantities together with default
  units for them. A Unit System has

  (required) handle bunits[7]; default units for 7 base
    quantities.

  (optional) char *name; name of this unit system

  (optional) handle dunits*; default units for any derived
    quantities that are NOT deriveable from bunits.

  Examples:

  International System of Units (SI) is an example of a system of
    units where default units for base quantities are meters,
    kilograms, seconds, amperes, kelvin, mole, candela. From
    these, derived units for velocity, for example,
    meters/second, current density, amperes/meter^2, etc.

  A and B Divisions at LLNL use default units such as 10
    nanoseconds for time, centimeters for length, etc.

    Notes: If in a given iMesh/iField instance we ever wind up
    mixing field objects created by different applications, then
    it is very likely going to be necessary to be able to
    explicitly set units for any given field (e.g. use something
    other than the default determined by some Unit System).
    This will be especially true if in any given iMesh/iField
    instance we support only one 'active' Unit System.

  It might be nice to allow applications to refer to systems of
  units by name rather than having to create them from scratch as
  needed.

\section Other Notes and Comments

   MCM: 'Precision' is potentially confusing here, particularly since this
   list is not either all floating point formats or all integer formats. Shouldn't
   this list really come out of iBase anyways and/or shouldn't we extend iBase
   to include types needed by iField?

   MCM: The notion expressed here seems to be a kind of hybrid of the data type
   as well as operations it supports. None of these are really algebraic types
   in the sense that I introduced that term for this API. 

   Though I see the possible value in distinguishing between dofs obeying
   boolean algebra semantics and those obeying say complex number semantics,
   that isn't what I had originally introduced the notion of algebraic type
   to support.

   At the same time, we're missing entirely the notion of vector valued and
   tensor valued fields here as well as 'scalar' which all the examples below
   are really instances of certain classes of algebraic types of 'scalar'.

   If it was just a matter of distinguishing between scalar, vector and tensor
   type things, we'd simply declare all fields to be 'tensor valued' and just
   specify the tensor order (which we do anyways). But, that then misses the
   nuances of integral valued vs. real valued vs. complex valued (what about
   integral-complex? digital signal processing codes have to deal with this
   issue). 

   Also, for something like 'complex', that involves 2 things, not one. And,
   those things could be real-imaginary pairs (cartesian) or magnitude-phase
   pairs (polar). In other words, there is a tiny little 'coordinate system'
   that defines how values of the field (thats the range set of the domain-range
   pair of sets needed to define a function) are interpreted. For a quaternion,
   there would be 4 things. I mention quaternion only to help in understanding
   how to properly generalize and abstract this. Even for a single thing (number),
   there is a difference between just that number and that number multiplying
   a vector, say the unit X direction vector. One is a scalar, the other is a
   vector even though both involve only a single number, the magnitude in the X
   direction. In other words, 'a' is different from 'a multiplying x-hat'. The
   latter is what you'd call the 'component' of some vector field and that is
   qualitatively different from something that is truly a scalar value'd field.

   The issue with whether or not the values of a field take on only discrete
   (integral or boolean even) values is a statement about the range set
   (of the domain<->range pair that we're defining a function between. It says
   whether or not the set is 'continuous' or 'discrete' as in the difference
   between the real number line and the integers.

   I think this best way to handle these notions is a slightly more complex API
   for defining algebraic types similar to the DFK api already being developed.

   CFOG reply:
   Mark, I don't see (at least yet) why the combination of tensor order,
   precision, and algebraic type doesn't cover the cases you discuss.  I
   agree that there may be some overlap.  Would be be better off keeping
   precision (bool, int, float, double, quad, ???) and having this
   variable tell us only real/complex/??? ?  I think I'd prefer not to
   roll the two together (nearly doubling the number of precisions).  In
   any event, all three of the pieces of info in my first sentence are
   going to be needed to decide on storage. 

   Before we go making this a user-extensible kind of thing, I'd
   definitely want to see real, concrete use cases that aren't covered
   by something simple.

\section Define possible storage hints for a tensor
  
   Depending on application context, the optimal way to store field data
   may differ.  For instance, a flow solver may prefer dofs for
   each entity to be stored contiguously (in an array of all dofs) for
   efficient access while an adaptation scheme may prefer dofs to be
   associated directly with entities so that dofs can be created and
   destroyed with entities.
  
   The application can provide suggestions to an implementation about
   what dof storage pattern is likely to be advantageous using the
   values enumerated here.  An implementation may or may not choose to
   take advantage of the hint.  Either way, the implementation is still
   responsible for correctly responding to all iField queries.
  
   Values of the enumeration are currently iField_BLOCKED,
   iField_INTERLEAVED, iField_MIXED, and iField_PER_ENTITY. 
  
   Block implies that a scalar temperature (T) and velocity (u,v,w) in a
   compound field are desired to be stored as: 
      (TTTTTT...uuuuuu...vvvvvv...wwwwww...).
   For a simple velocity field, omit the T's in this example.
  
   Interleaved implies that a scalar temperature (T) and velocity
   (u,v,w) in a compound field are desired to be stored as: 
      (TuvwTuvwTuvwTuvwTuvwTuvw...).
   For a simple velocity field, omit the T's in this example.
  
   Mixed implies that a scalar temperature (T) and velocity (u,v,w) in a
   compound field are desired to be stored as: 
       (TTTTTT...uvwuvwuvwuvwuvwuvw...).
   For simple fields, this option is not allowed.
  
   Per entity implies that the app wants data associated with entities
   to support mesh modification / adaptation.  (Yes, these can be
   supported using array-based storage, but it's harder that way.)
  
   MCM or FD: We can't support the case where a caller has 3 separate
   arrays one for X, one for Y and one for Z. That seems like way too
   common a case to NOT handle. In fact, for many rank=1 or 2 tensor
   fields, this is a real common storage paradigm. All we really need to
   do to handle this case is add iField_SEPARATE.
  
   CFOG reply:  For data stored internally in the iField implementation,
   how does iField_SEPARATE differ functionally from iField_BLOCKED?
   Whether I use one array or separate arrays, I'm pretty sure that
   minor pointer magic in the implementation will result in identical
   machine code.

   CFOG 9/26/10  I can definitely see wanting to iterate separately over
   entities with different dfuncs.  Should these be separate calls to
   create iterators, or additional args to the existing iterator calls?

   Also, IMO, apps are responsible for ensuring sensible behavior in the
   transitional cases (entities with some faces p-refined but not all,
   say.  If nothing else, how this is dealt with depends on what dfuncs
   you combine, so the implementation needs guidance anyway.  Caveat
   user.

   NOTES / QUESTIONS / TOPICS FOR DISCUSSION:
  
   Placeholder: At a minimum, we need to be able to represent tensor
   product spaces $(\Re^{n}\times\Re)$ --- space-time is the most common
   but by no means only example.  Regardless of details, this requires
   more thought.
  
   Placeholder: Variable order (p-refinement)
  
   Placeholder: Thread safety: atomic += and lock/unlock?  Or strictly
   up to implementation?
  
   Placeholder:  Returning mixed real / int data.  Up-convert to double?
  
   Placeholder:  dfunc kernel support:  iBase_Type, plus indication of
   dual? 
  
   1.  Being able to express relationships between fields is a useful
       capability, even if iField does nothing more than store the
       relationship and emit it on demand.  We discussed differential
       and algebraic operators; others are possible as well.  Evaluation
       of these operators is something that is probably best done as a
       service (implementations could of course choose to implement this
       service directly, making operations on their own fields more
       efficient...).  We tentatively concluded that this capability is
       unlikely to affect the core of the iField API and may even be
       orthogonal to most or all of the rest of the API.
  
   4.  Coordinate consistency: I am thinking there is an argument for
       creating a (non-coordinate) field using a mesh topology and a
       coordinate field.  The latter would carry info about cartesian
       vs. polar, etc.  Physical coordinates passed to that field (and
       its component tensors and dfuncs, where appropriate) would be in
       that coordinate systems.  There are some obvious advantages to
       enforcing coordinate system consistency in this way, but also
       some possible disadvantages.  Any thoughts on this?
  
   5.  We may want functions to be able to iterate over all dof clusters
       (dofs associated with the same entity) for things like changes of
       variables and/or coordinate system. Otherwise, you either have to
       retrieve all the dofs at once to do this, or duplicate work by
       doing the conversion (in a copy mode) multiple times for shared
       dofs, or keep track of which dofs you have done.
 ******************************************************************************/

#ifdef __cplusplus
} /*  extern "C"  */
#endif

#endif /* defined(iField_H) */
