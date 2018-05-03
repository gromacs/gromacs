/* ************************************************************************
 * Copyright 2013 Advanced Micro Devices, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ************************************************************************/


/*! @file clAmdFft.h
 * /note clAmdFft.h is a deprecated header file.  
 * This header is provided to help projects that were written with the older clAmdFft codebase, to help them 
 * port to the new API at their own schedule.  It will not be maintained or updated, and will be removed after 
 * a reasonable amount of time has passed.  All new code should be written against clFFT.h.  
 * Older projects should migrate to the new header at their earliest convenience.
 */

#pragma once
#if !defined( CLAMDFFT_DOTH )
#define CLAMDFFT_DOTH

#include "clFFT.h"

/* The following header defines a fixed version number as this header is deprecated and won't be updated */
#include "clAmdFft.version.h"

/*	In general, you can not use namespaces for strict C compliance, so we prefix our public accessible names
 *	with the string clAmdFft
 */

/*	All functions will return pre-defined error codes, and will NOT throw exceptions to the caller
 */

/*!  @brief clAmdFft error codes definition, incorporating OpenCL error definitions
 *
 *   This enumeration is a superset of the OpenCL error codes.  For example, CL_OUT_OF_HOST_MEMORY,
 *   which is defined in cl.h is aliased as CLFFT_OUT_OF_HOST_MEMORY.  The set of basic OpenCL
 *   error codes is extended to add extra values specific to the clAmdFft package.
 */
typedef enum clfftStatus_ clAmdFftStatus;

/*!  @brief The dimension of the input and output buffers that will be fed into all FFT transforms */
typedef enum clfftDim_ clAmdFftDim;

/*!  @brief These are the expected layouts of the buffers */
typedef enum clfftLayout_ clAmdFftLayout;

/*!  @brief This is the expected precision of each FFT.
 */
typedef enum clfftPrecision_ clAmdFftPrecision;

/*!  @brief What is the expected direction of each FFT, time or the frequency domains */
typedef enum clfftDirection_ clAmdFftDirection;

/*!  @brief Are the input buffers overwritten with the results */
typedef enum clfftResultLocation_ clAmdFftResultLocation;

/*! @brief This determines whether the result is returned in original order. It is valid only for
dimensions greater than 1. */
typedef enum clfftResultTransposed_ clAmdFftResultTransposed;

/*! @brief Data structure that can be passed to clAmdFftSetup() to control the behavior of the FFT runtime
 *  @details This structure contains values that can be initialized before instantiation of the FFT runtime
 *  with ::clAmdFftSetup().  To initialize this structure, pass a pointer to a user struct to ::clAmdFftInitSetupData( ),
 *  which will clear the structure and set the version member variables to the current values.
 */
typedef struct clfftSetupData_ clAmdFftSetupData;

/*!  @brief An abstract handle to the object that represents the state of the FFT(s) */
typedef clfftPlanHandle clAmdFftPlanHandle;

#ifdef __cplusplus
extern "C" {
#endif
	/*! @brief Initialize an clAmdFftSetupData struct for the client
	 *  @details clAmdFftSetupData is passed to clAmdFftSetup to control behavior of the FFT runtime
	 *  @param[out] setupData Data structure is cleared, initialized with version information and default values
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftInitSetupData( clAmdFftSetupData* setupData )
	{
		return clfftInitSetupData( setupData );
	}

	/*! @brief Initialize internal FFT resources.
	 *  @details AMD's FFT implementation caches kernels, programs and buffers for its internal use.
	 *  @param[in] setupData Data structure that can be passed into the setup routine to control FFT generation behavior
	 * 	and debug functionality
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftSetup( const clAmdFftSetupData* setupData )
	{
		return clfftSetup( setupData );
	}

	/*! @brief Release all internal resources.
	 *  @details Call when client is done with this FFT library, allowing the library to destroy all resources it has cached
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftTeardown( )
	{
		return clfftTeardown( );
	}

	/*! @brief Query the FFT library for version information
	 *  @details Return the major, minor and patch version numbers associated with this FFT library
	 *  @param[out] major Major functionality change
	 *  @param[out] minor Minor functionality change
	 *  @param[out] patch Bug fixes, documentation changes, no new features introduced
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftGetVersion( cl_uint* major, cl_uint* minor, cl_uint* patch )
	{
		return clfftGetVersion( major, minor, patch );
	}

	/*! @brief Create a plan object initialized entirely with default values.
	 *  @details A plan is a repository of state for calculating FFT's.  Allows the runtime to pre-calculate kernels, programs
	 * 	and buffers and associate them with buffers of specified dimensions.
	 *  @param[out] plHandle Handle to the newly created plan
	 *  @param[in] context Client is responsible for providing an OpenCL context for the plan
	 *  @param[in] dim The dimensionality of the FFT transform; describes how many elements are in the array
	 *  @param[in] clLengths An array of lengths, of size 'dim'.  Each value describes the length of additional dimensions
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftCreateDefaultPlan( clAmdFftPlanHandle* plHandle, cl_context context, const clAmdFftDim dim,
								const size_t* clLengths )
	{
		return clfftCreateDefaultPlan( plHandle, context, dim, clLengths );
	}

	/*! @brief Create a copy of an existing plan.
	 *  @details This API allows a client to create a new plan based upon an existing plan.  This is a convenience function
	 *  provided for quickly creating plans that are similar, but may differ slightly.
	 *  @param[out] out_plHandle Handle to the newly created plan that is based on in_plHandle
	 *  @param[in] new_context Client is responsible for providing a new context for the new plan
	 *  @param[in] in_plHandle Handle to a plan to be copied, previously created
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftCopyPlan( clAmdFftPlanHandle* out_plHandle, cl_context new_context, clAmdFftPlanHandle in_plHandle )
	{
		return clfftCopyPlan( out_plHandle, new_context, in_plHandle );
	}

	/*! @brief Prepare the plan for execution.
	 *  @details After all plan parameters are set, the client has the option of 'baking' the plan, which tells the runtime that
	 *  no more changes to the plan's parameters are expected, and the OpenCL kernels should be compiled.  This optional function
	 *  allows the client application to perform this function when the application is being initialized instead of on the first
	 *  execution.
	 *  At this point, the clAmdFft runtime will apply all implimented optimizations, possibly including
	 *  running kernel experiments on the devices in the plan context.
	 *  <p>  Users should assume that this function will take a long time to execute.  If a plan is not baked before being executed,
	 *  users should assume that the first call to clAmdFftEnqueueTransform will take a long time to execute.
	 *  <p>  If any significant parameter of a plan is changed after the plan is baked (by a subsequent call to one of
	 *  the clAmdFftSetPlan____ functions), that will not be considered an error.  Instead, the plan will revert back to
	 *  the unbaked state, discarding the benefits of the baking operation.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] numQueues Number of command queues in commQueueFFT; 0 is a valid value, in which case client does not want
	 * 	the runtime to run load experiments and only pre-calculate state information
	 *  @param[in] commQueueFFT An array of cl_command_queues created by the client; the command queues must be a proper subset of
	 * 	the devices included in the plan context
	 *  @param[in] pfn_notify A function pointer to a notification routine. The notification routine is a callback function that
	 *  an application can register and which will be called when the program executable has been built (successfully or unsuccessfully).
	 *  Currently, this parameter MUST be NULL or nullptr.
	 *  @param[in] user_data Passed as an argument when pfn_notify is called.
	 *  Currently, this parameter MUST be NULL or nullptr.
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftBakePlan( clAmdFftPlanHandle plHandle, cl_uint numQueues, cl_command_queue* commQueueFFT,
							void (CL_CALLBACK *pfn_notify)(clAmdFftPlanHandle plHandle, void *user_data), void* user_data )
	{
		return clfftBakePlan( plHandle, numQueues, commQueueFFT, pfn_notify, user_data );
	}

	/*! @brief Release the resources of a plan.
	 *  @details A plan may include kernels, programs and buffers associated with it that consume memory.  When a plan
	 *  is not needed anymore, the client should release the plan.
	 *  @param[in,out] plHandle Handle to a plan previously created
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftDestroyPlan( clAmdFftPlanHandle* plHandle )
	{
		return clfftDestroyPlan( plHandle );
	}

	/*! @brief Retrieve the OpenCL context of a previously created plan.
	 *  @details User should pass a reference to an cl_context variable, which will be changed to point to a
	 *  context set in the specified plan.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] context Reference to user allocated cl_context, which will point to context set in plan
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftGetPlanContext( const clAmdFftPlanHandle plHandle, cl_context* context )
	{
		return clfftGetPlanContext( plHandle, context );
	}

	/*! @brief Retrieve the floating point precision of the FFT data
	 *  @details User should pass a reference to an clAmdFftPrecision variable, which will be set to the
	 *  precision of the FFT complex data in the plan.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] precision Reference to user clAmdFftPrecision enum
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftGetPlanPrecision( const clAmdFftPlanHandle plHandle, clAmdFftPrecision* precision )
	{
		return clfftGetPlanPrecision( plHandle, precision );
	}

	/*! @brief Set the floating point precision of the FFT data
	 *  @details Set the plan property which will be the precision of the FFT complex data in the plan.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] precision Reference to user clAmdFftPrecision enum
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftSetPlanPrecision( clAmdFftPlanHandle plHandle, clAmdFftPrecision precision )
	{
		return clfftSetPlanPrecision( plHandle, precision );
	}

	/*! @brief Retrieve the scaling factor that should be applied to the FFT data
	 *  @details User should pass a reference to an cl_float variable, which will be set to the
	 *  floating point scaling factor that will be multiplied across the FFT data.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dir Which direction does the scaling factor apply to
	 *  @param[out] scale Reference to user cl_float variable
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftGetPlanScale( const clAmdFftPlanHandle plHandle, clAmdFftDirection dir, cl_float* scale )
	{
		return clfftGetPlanScale( plHandle, dir, scale );
	}

	/*! @brief Set the scaling factor that should be applied to the FFT data
	 *  @details Set the plan property which will be the floating point scaling factor that will be
	 *  multiplied across the FFT data.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dir Which direction does the scaling factor apply to
	 *  @param[in] scale Reference to user cl_float variable
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftSetPlanScale( clAmdFftPlanHandle plHandle, clAmdFftDirection dir, cl_float scale )
	{
		return clfftSetPlanScale( plHandle, dir, scale );
	}

	/*! @brief Retrieve the number of discrete arrays that this plan can handle concurrently
	 *  @details User should pass a reference to an cl_uint variable, which will be set to the
	 *  number of discrete arrays (1D or 2D) that will be batched together for this plan
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] batchSize How many discrete number of FFT's are to be performed
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftGetPlanBatchSize( const clAmdFftPlanHandle plHandle, size_t* batchSize )
	{
		return clfftGetPlanBatchSize( plHandle, batchSize );
	}

	/*! @brief Set the number of discrete arrays that this plan can handle concurrently
	 *  @details Set the plan property which will be set to the number of discrete arrays (1D or 2D)
	 *  that will be batched together for this plan
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] batchSize How many discrete number of FFT's are to be performed
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftSetPlanBatchSize( clAmdFftPlanHandle plHandle, size_t batchSize )
	{
		return clfftSetPlanBatchSize( plHandle, batchSize );
	}

	/*! @brief Retrieve the dimensionality of FFT's to be transformed in the plan
	 *  @details Queries a plan object and retrieves the dimensionality that the plan is set for.  A size is returned to
	 *  help the client allocate the proper storage to hold the dimensions in a further call to clAmdFftGetPlanLength
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] dim The dimensionality of the FFT's to be transformed
	 *  @param[out] size Value used to allocate an array to hold the FFT dimensions.
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftGetPlanDim( const clAmdFftPlanHandle plHandle, clAmdFftDim* dim, cl_uint* size )
	{
		return clfftGetPlanDim( plHandle, dim, size );
	}

	/*! @brief Set the dimensionality of FFT's to be transformed by the plan
	 *  @details Set the dimensionality of FFT's to be transformed by the plan
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dim The dimensionality of the FFT's to be transformed
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftSetPlanDim( clAmdFftPlanHandle plHandle, const clAmdFftDim dim )
	{
		return clfftSetPlanDim( plHandle, dim );
	}

	/*! @brief Retrieve the length of each dimension of the FFT
	 *  @details User should pass a reference to a size_t array, which will be set to the
	 *  length of each discrete dimension of the FFT
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dim The dimension of the length parameters; describes how many elements are in the array
	 *  @param[out] clLengths An array of lengths, of size 'dim'.  Each array value describes the length of each dimension
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftGetPlanLength( const clAmdFftPlanHandle plHandle, const clAmdFftDim dim, size_t* clLengths )
	{
		return clfftGetPlanLength( plHandle, dim, clLengths );
	}

	/*! @brief Set the length of each dimension of the FFT
	 *  @details Set the plan property which will be the length of each discrete dimension of the FFT
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dim The dimension of the length parameters; describes how many elements are in the array
	 *  @param[in] clLengths An array of lengths, of size 'dim'.  Each value describes the length of additional dimensions
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftSetPlanLength( clAmdFftPlanHandle plHandle, const clAmdFftDim dim, const size_t* clLengths )
	{
		return clfftSetPlanLength( plHandle, dim, clLengths );
	}

	/*! @brief Retrieve the distance between consecutive elements for input buffers in a dimension.
	 *  @details Depending on how the dimension is set in the plan (for 2D or 3D FFT's), strideY or strideZ can be safely
	 *  ignored
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dim The dimension of the stride parameters; describes how many elements are in the array
	 *  @param[out] clStrides An array of strides, of size 'dim'.
	 */
	__inline clAmdFftStatus clAmdFftGetPlanInStride( const clAmdFftPlanHandle plHandle, const clAmdFftDim dim, size_t* clStrides )
	{
		return clfftGetPlanInStride( plHandle, dim, clStrides );
	}

	/*! @brief Set the distance between consecutive elements for input buffers in a dimension.
	 *  @details Set the plan properties which will be the distance between elements in a given dimension
	 *  (units are in terms of clAmdFftPrecision)
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dim The dimension of the stride parameters; describes how many elements are in the array
	 *  @param[in] clStrides An array of strides, of size 'dim'. Usually strideX=1 so that successive elements in the first dimension are stored contiguously.
	 * 	Typically strideY=LenX, strideZ=LenX*LenY such that successive elements in the second and third dimensions are stored in packed format.
	 *  See  @ref DistanceStridesandPitches for details.
	 */
	__inline clAmdFftStatus clAmdFftSetPlanInStride( clAmdFftPlanHandle plHandle, const clAmdFftDim dim, size_t* clStrides )
	{
		return clfftSetPlanInStride( plHandle, dim, clStrides );
	}

	/*! @brief Retrieve the distance between consecutive elements for output buffers in a dimension.
	 *  @details Depending on how the dimension is set in the plan (for 2D or 3D FFT's), strideY or strideZ can be safely
	 *  ignored
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dim The dimension of the stride parameters; describes how many elements are in the array
	 *  @param[out] clStrides An array of strides, of size 'dim'.
	 */
	__inline clAmdFftStatus clAmdFftGetPlanOutStride( const clAmdFftPlanHandle plHandle, const clAmdFftDim dim, size_t* clStrides )
	{
		return clfftGetPlanOutStride( plHandle, dim, clStrides );
	}

	/*! @brief Set the distance between consecutive elements for output buffers in a dimension.
	 *  @details Set the plan properties which will be the distance between elements in a given dimension
	 *  (units are in terms of clAmdFftPrecision)
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dim The dimension of the stride parameters; describes how many elements are in the array
	 *  @param[in] clStrides An array of strides, of size 'dim'.  Usually strideX=1 so that successive elements in the first dimension are stored contiguously.
	 * 	Typically strideY=LenX, strideZ=LenX*LenY such that successive elements in the second and third dimensions are stored in packed format.
	 *  @sa clAmdFftSetPlanInStride
	 */
	__inline clAmdFftStatus clAmdFftSetPlanOutStride( clAmdFftPlanHandle plHandle, const clAmdFftDim dim, size_t* clStrides )
	{
		return clfftSetPlanOutStride( plHandle, dim, clStrides );
	}

	/*! @brief Retrieve the distance between Array objects
	 *  @details Pitch is the distance between each discrete array object in an FFT array. This is only used
	 *  for 'array' dimensions in clAmdFftDim; see clAmdFftSetPlanDimension (units are in terms of clAmdFftPrecision)
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] iDist The distance between the beginning elements of the discrete array objects in memory on input.
	 *  For contiguous arrays in memory, iDist=(strideX*strideY*strideZ)
	 *  @param[out] oDist The distance between the beginning elements of the discrete array objects in memory on output.
	 *  For contiguous arrays in memory, oDist=(strideX*strideY*strideZ)
	 */
	__inline clAmdFftStatus clAmdFftGetPlanDistance( const clAmdFftPlanHandle plHandle, size_t* iDist, size_t* oDist )
	{
		return clfftGetPlanDistance( plHandle, iDist, oDist );
	}

	/*! @brief Set the distance between Array objects
	 *  @details Pitch is the distance between each discrete array object in an FFT array. This is only used
	 *  for 'array' dimensions in clAmdFftDim; see clAmdFftSetPlanDimension (units are in terms of clAmdFftPrecision)
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] iDist The distance between the beginning elements of the discrete array objects in memory on input.
	 *  For contiguous arrays in memory, iDist=(strideX*strideY*strideZ)
	 *  @param[out] oDist The distance between the beginning elements of the discrete array objects in memory on output.
	 *  For contiguous arrays in memory, oDist=(strideX*strideY*strideZ)
	 */
	__inline clAmdFftStatus clAmdFftSetPlanDistance( clAmdFftPlanHandle plHandle, size_t iDist, size_t oDist )
	{
		return clfftSetPlanDistance( plHandle, iDist, oDist );
	}

	/*! @brief Retrieve the expected layout of the input and output buffers
	 *  @details Output buffers can be filled with either hermitian or complex numbers.  Complex numbers can be stored
	 *  in various layouts; this informs the FFT engine what layout to produce on output
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] iLayout Indicates how the input buffers are laid out in memory
	 *  @param[out] oLayout Indicates how the output buffers are laid out in memory
	 */
	__inline clAmdFftStatus clAmdFftGetLayout( const clAmdFftPlanHandle plHandle, clAmdFftLayout* iLayout, clAmdFftLayout* oLayout )
	{
		return clfftGetLayout( plHandle, iLayout, oLayout );
	}

	/*! @brief Set the expected layout of the input and output buffers
	 *  @details Output buffers can be filled with either hermitian or complex numbers.  Complex numbers can be stored
	 *  in various layouts; this informs the FFT engine what layout to produce on output
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] iLayout Indicates how the input buffers are laid out in memory
	 *  @param[in] oLayout Indicates how the output buffers are laid out in memory
	 */
	__inline clAmdFftStatus clAmdFftSetLayout( clAmdFftPlanHandle plHandle, clAmdFftLayout iLayout, clAmdFftLayout oLayout )
	{
		return clfftSetLayout( plHandle, iLayout, oLayout );
	}

	/*! @brief Retrieve whether the input buffers are going to be overwritten with results
	 *  @details If the setting is to do an in-place transform, the input buffers are overwritten with the results of the
	 *  transform.  If the setting is for out-of-place transforms, the engine knows to look for separate output buffers
	 *  on the Enqueue call.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] placeness Tells the FFT engine to clobber the input buffers or to expect output buffers for results
	 */
	__inline clAmdFftStatus clAmdFftGetResultLocation( const clAmdFftPlanHandle plHandle, clAmdFftResultLocation* placeness )
	{
		return clfftGetResultLocation( plHandle, placeness );
	}

	/*! @brief Set whether the input buffers are going to be overwritten with results
	 *  @details If the setting is to do an in-place transform, the input buffers are overwritten with the results of the
	 *  transform.  If the setting is for out-of-place transforms, the engine knows to look for separate output buffers
	 *  on the Enqueue call.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] placeness Tells the FFT engine to clobber the input buffers or to expect output buffers for results
	 */
	__inline clAmdFftStatus clAmdFftSetResultLocation( clAmdFftPlanHandle plHandle, clAmdFftResultLocation placeness )
	{
		return clfftSetResultLocation( plHandle, placeness );
	}

	/*! @brief Retrieve the final transpose setting of a muti-dimensional FFT
	 *  @details A multi-dimensional FFT typically transposes the data several times during calculation.  If the client
	 *  does not care about the final transpose to put data back in proper dimension, the final transpose can be skipped
	 *  for possible speed improvements
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] transposed Parameter specifies whether the final transpose can be skipped
	 */
	__inline clAmdFftStatus clAmdFftGetPlanTransposeResult( const clAmdFftPlanHandle plHandle, clAmdFftResultTransposed * transposed )
	{
		return clfftGetPlanTransposeResult( plHandle, transposed );
	}

	/*! @brief Set the final transpose setting of a muti-dimensional FFT
	 *  @details A multi-dimensional FFT typically transposes the data several times during calculation.  If the client
	 *  does not care about the final transpose to put data back in proper dimension, the final transpose can be skipped
	 *  for possible speed improvements
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] transposed Parameter specifies whether the final transpose can be skipped
	 */
	__inline clAmdFftStatus clAmdFftSetPlanTransposeResult( clAmdFftPlanHandle plHandle, clAmdFftResultTransposed transposed )
	{
		return clfftSetPlanTransposeResult( plHandle, transposed );
	}

	/*! @brief Get buffer size (in bytes), which may be needed internally for an intermediate buffer
	 *  @details Very large FFT transforms may need multiple passes, and the operation would need a temporary buffer to hold
	 *  intermediate results. This function is only valid after the plan is baked, otherwise an invalid operation error
	 *  is returned. If buffersize returns as 0, the runtime needs no temporary buffer.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[out] buffersize Size in bytes for intermediate buffer
	 */
	__inline clAmdFftStatus clAmdFftGetTmpBufSize( const clAmdFftPlanHandle plHandle, size_t* buffersize )
	{
		return clfftGetTmpBufSize( plHandle, buffersize );
	}

	/*! @brief Enqueue an FFT transform operation, and return immediately (non-blocking)
	 *  @details This transform API is the function that actually computes the FFT transfrom. It is non-blocking as it
	 *  only enqueues the OpenCL kernels for execution. The synchronization step has to be managed by the user.
	 *  @param[in] plHandle Handle to a plan previously created
	 *  @param[in] dir Forwards or backwards transform
	 *  @param[in] numQueuesAndEvents Number of command queues in commQueues; number of expected events to be returned in outEvents
	 *  @param[in] commQueues An array of cl_command_queues created by the client; the command queues must be a proper subset of
	 * 	the devices included in the plan context
	 *  @param[in] numWaitEvents Specify the number of elements in the eventWaitList array
	 *  @param[in] waitEvents Events that this transform should wait to complete before executing on the device
	 *  @param[out] outEvents The runtime fills this array with events corresponding 1 to 1 with the input command queues passed
	 *	in commQueues.  This parameter can be NULL or nullptr, in which case client is not interested in receiving notifications
	 *	when transforms are finished, otherwise if not NULL the client is responsible for allocating this array, with at least
	 *	as many elements as specified in numQueuesAndEvents.
	 *  @param[in] inputBuffers An array of cl_mem objects that contain data for processing by the FFT runtime.  If the transform
	 *  is in place, the FFT results will overwrite the input buffers
	 *  @param[out] outputBuffers An array of cl_mem objects that will store the results of out of place transforms.  If the transform
	 *  is in place, this parameter may be NULL or nullptr.  It is completely ignored
	 *  @param[in] tmpBuffer A cl_mem object that is reserved as a temporary buffer for FFT processing. If clTmpBuffers is NULL or nullptr,
	 *  and the runtime needs temporary storage, an internal temporary buffer will be created on the fly managed by the runtime.
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	__inline clAmdFftStatus clAmdFftEnqueueTransform(
												clAmdFftPlanHandle plHandle,
												clAmdFftDirection dir,
												cl_uint numQueuesAndEvents,
												cl_command_queue* commQueues,
												cl_uint numWaitEvents,
												const cl_event* waitEvents,
												cl_event* outEvents,
												cl_mem* inputBuffers,
												cl_mem* outputBuffers,
												cl_mem tmpBuffer
												)
	{
		return clfftEnqueueTransform( plHandle, dir, numQueuesAndEvents, commQueues, numWaitEvents, waitEvents, outEvents, 
			inputBuffers, outputBuffers, tmpBuffer );
	}

#ifdef __cplusplus
}
#endif

#endif
