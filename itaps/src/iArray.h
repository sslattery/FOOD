#ifndef _ITAPS_iArray
#define _ITAPS_iArray

#include "iBase.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Enumeration of ordering of array dimensions.
 */
enum iArray_DimensionOrder {
    iArray_DimensionOrder_MIN = 0,
    iArray_C_ORDER = iArray_DimensionOrder_MIN,
    iArray_FORTRAN_ORDER,
    iArray_DimensionOrder_MAX = iArray_FORTRAN_ORDER
};

/*!
 * \brief Describe the status of the array.
 */
void iArray_getDescription( iBase_EntityHandle array,
			    char* descr,
			    int descr_len );

/*!
 * \brief Allocate a new array of arbitray data type.
 */
void iArray_allocateArray( const char* options,
			   iBase_EntityHandle* array,
			   const int storage_order,
			   const int dimension_order,
			   const int* dimensions,
			   const int num_dimension,
			   const int* data_allocated,
			   int* err,
			   int options_len );

/*!
 * \brief Create an array of arbitrary data type using existing data.
 */
void iArray_createArray(  const char* options,
			  iBase_EntityHandle* array,
			  const int storage_order,
			  const int dimension_order,
			  const int* dimensions,
			  const int num_dimension,
			  void* data,
			  int* data_allocated,
			  int* data_size,
			  const int copy,
			  int* err,
			  int options_len );

/*! 
 * \brief Destroy an array.
 */
void iArray_dtor( iBase_EntityHandle array,
		  int* error );

/*!
 * \brief Get the size of an array in bytes.
 */
void iArray_getSize( const iBase_EntityHandle array,
		     int* size,
		     int* err );

/*!
 * \brief Get the number of dimensions of an array.
 */
void iArray_getNumDimensions( const iBase_EntityHandle array,
			      int* num_dimensions,
			      int *err );

/*!
 * \brief Get the size of the dimensions of an array.
 */
void iArray_getDimensions( const iBase_EntityHandle array,
			   int* dimensions,
			   int* dimensions_allocated,
			   int* dimensions_size,
			   int* err );

/*!
 * \brief Get the size of a particular dimension of an array.
 */
void iArray_getDimension( const iBase_EntityHandle array,
			  const int dimension,
			  int* dimension_size,
			  int* error );
			   
/*!
 * \brief Set the data in an array for an arbitrary type.
 */
void iArray_setData( iBase_EntityHandle array,
		     const int storage_order,
		     const void* data,
		     const int data_size,
		     const int copy,
		     int* err );

/*!
 * \brief Get the data in an array for an arbitrary type.
 */
void iArray_getData( const iBase_EntityHandle array,
		     const int storage_order,
		     void* data,
		     int* data_allocated,
		     int* data_size,
		     int* err );

/*!
 * \brief Set an individual value in the array.
 */
void iArray_setValue( iBase_EntityHandle array,
		      const int* indices,
		      const int num_indices,
		      const void* value,
		      const int* value_size,
		      int* err );

/*!
 * \brief Get an individual value in the array
 */
void iArray_getValue( const iBase_EntityHandle array,
		      const int* indices,
		      const int num_indices,
		      void* value,
		      int* value_size,
		      int* err );

/*!
 * \brief Set a slice of data in the array.
 */
void iArray_setSlice( iBase_EntityHandle array,
		      const int storage_order,
		      const int* slice_dimensions,
		      const int* slice_dimension_indices,
		      const int num_slice_parameters,
		      const void* data,
		      const int data_size,
		      int* err );

/*!
 * \brief Get a slice of data from the array.
 */
void iArray_getSlice( const iBase_EntityHandle array,
		      const int storage_order,
		      const int* slice_dimensions,
		      const int* slice_dimension_indices,
		      const int num_slice_parameters,
		      void* data,
		      int* data_allocated,
		      int* data_size,
		      int* err );

#ifdef __cplusplus
} // end extern "C"
#endif

#endif // end _ITAPS_iArray

