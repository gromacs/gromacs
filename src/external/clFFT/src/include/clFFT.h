/* ************************************************************************
 * Copyright 2013-2015 Advanced Micro Devices, Inc.
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


/*! @file clFFT.h
 * clFFT.h defines all the public interfaces and types that are used by clFFT clients
 * This is the only public header file that should be consumed by clFFT clients.  It is written to adhere to native "C"
 * interfaces to make clFFT library as portable as possible; it should be callable from C, C++, .NET and Fortran,
 * either with the proper linking or using wrapper classes.
 *
 */

#pragma once
#if !defined( CLFFT_H )
#define CLFFT_H

#if defined(__APPLE__) || defined(__MACOSX)
	#include <OpenCL/cl.h>
#else
	#include <CL/cl.h>
#endif

#include "clFFT.version.h"

/*! This preprocessor definition is the standard way to export APIs
 *  from a DLL simpler. All files within this DLL are compiled with the CLFFT_EXPORTS
 *  symbol defined on the command line. This symbol must not be defined on any project
 *  that uses this DLL. This ensures source files of any other project that include this file see
 *  clfft functions as being imported from a DLL, whereas the DLL sees symbols
 *  defined with this macro as being exported.
 */
#if defined( _WIN32 )
	#if !defined( __cplusplus )
		#define inline __inline
	#endif

    #if defined( CLFFT_STATIC )
        #define CLFFTAPI
    #elif defined( CLFFT_EXPORTS )
        #define CLFFTAPI __declspec( dllexport )
    #else
        #define CLFFTAPI __declspec( dllimport )
    #endif
#else
	#define CLFFTAPI
#endif

/*	In general, you cannot use namespaces for strict C compliance, so we prefix our public accessible names
 *	with the string clfft
 */

/*	All functions return pre-defined error codes, and do NOT throw exceptions to the caller.
 */

/*!  @brief clfft error codes definition(incorporating OpenCL error definitions)
 *
 *   This enumeration is a superset of the OpenCL error codes.  For example, CL_OUT_OF_HOST_MEMORY,
 *   which is defined in cl.h is aliased as CLFFT_OUT_OF_HOST_MEMORY.  The set of basic OpenCL
 *   error codes is extended to add extra values specific to the clfft package.
 */
enum clfftStatus_
{
	CLFFT_INVALID_GLOBAL_WORK_SIZE			= CL_INVALID_GLOBAL_WORK_SIZE,
	CLFFT_INVALID_MIP_LEVEL					= CL_INVALID_MIP_LEVEL,
	CLFFT_INVALID_BUFFER_SIZE				= CL_INVALID_BUFFER_SIZE,
	CLFFT_INVALID_GL_OBJECT					= CL_INVALID_GL_OBJECT,
	CLFFT_INVALID_OPERATION					= CL_INVALID_OPERATION,
	CLFFT_INVALID_EVENT						= CL_INVALID_EVENT,
	CLFFT_INVALID_EVENT_WAIT_LIST			= CL_INVALID_EVENT_WAIT_LIST,
	CLFFT_INVALID_GLOBAL_OFFSET				= CL_INVALID_GLOBAL_OFFSET,
	CLFFT_INVALID_WORK_ITEM_SIZE			= CL_INVALID_WORK_ITEM_SIZE,
	CLFFT_INVALID_WORK_GROUP_SIZE			= CL_INVALID_WORK_GROUP_SIZE,
	CLFFT_INVALID_WORK_DIMENSION			= CL_INVALID_WORK_DIMENSION,
	CLFFT_INVALID_KERNEL_ARGS				= CL_INVALID_KERNEL_ARGS,
	CLFFT_INVALID_ARG_SIZE					= CL_INVALID_ARG_SIZE,
	CLFFT_INVALID_ARG_VALUE					= CL_INVALID_ARG_VALUE,
	CLFFT_INVALID_ARG_INDEX					= CL_INVALID_ARG_INDEX,
	CLFFT_INVALID_KERNEL					= CL_INVALID_KERNEL,
	CLFFT_INVALID_KERNEL_DEFINITION			= CL_INVALID_KERNEL_DEFINITION,
	CLFFT_INVALID_KERNEL_NAME				= CL_INVALID_KERNEL_NAME,
	CLFFT_INVALID_PROGRAM_EXECUTABLE		= CL_INVALID_PROGRAM_EXECUTABLE,
	CLFFT_INVALID_PROGRAM					= CL_INVALID_PROGRAM,
	CLFFT_INVALID_BUILD_OPTIONS				= CL_INVALID_BUILD_OPTIONS,
	CLFFT_INVALID_BINARY					= CL_INVALID_BINARY,
	CLFFT_INVALID_SAMPLER					= CL_INVALID_SAMPLER,
	CLFFT_INVALID_IMAGE_SIZE				= CL_INVALID_IMAGE_SIZE,
	CLFFT_INVALID_IMAGE_FORMAT_DESCRIPTOR	= CL_INVALID_IMAGE_FORMAT_DESCRIPTOR,
	CLFFT_INVALID_MEM_OBJECT				= CL_INVALID_MEM_OBJECT,
	CLFFT_INVALID_HOST_PTR					= CL_INVALID_HOST_PTR,
	CLFFT_INVALID_COMMAND_QUEUE				= CL_INVALID_COMMAND_QUEUE,
	CLFFT_INVALID_QUEUE_PROPERTIES			= CL_INVALID_QUEUE_PROPERTIES,
	CLFFT_INVALID_CONTEXT					= CL_INVALID_CONTEXT,
	CLFFT_INVALID_DEVICE					= CL_INVALID_DEVICE,
	CLFFT_INVALID_PLATFORM					= CL_INVALID_PLATFORM,
	CLFFT_INVALID_DEVICE_TYPE				= CL_INVALID_DEVICE_TYPE,
	CLFFT_INVALID_VALUE						= CL_INVALID_VALUE,
	CLFFT_MAP_FAILURE						= CL_MAP_FAILURE,
	CLFFT_BUILD_PROGRAM_FAILURE				= CL_BUILD_PROGRAM_FAILURE,
	CLFFT_IMAGE_FORMAT_NOT_SUPPORTED		= CL_IMAGE_FORMAT_NOT_SUPPORTED,
	CLFFT_IMAGE_FORMAT_MISMATCH				= CL_IMAGE_FORMAT_MISMATCH,
	CLFFT_MEM_COPY_OVERLAP					= CL_MEM_COPY_OVERLAP,
	CLFFT_PROFILING_INFO_NOT_AVAILABLE		= CL_PROFILING_INFO_NOT_AVAILABLE,
	CLFFT_OUT_OF_HOST_MEMORY				= CL_OUT_OF_HOST_MEMORY,
	CLFFT_OUT_OF_RESOURCES					= CL_OUT_OF_RESOURCES,
	CLFFT_MEM_OBJECT_ALLOCATION_FAILURE		= CL_MEM_OBJECT_ALLOCATION_FAILURE,
	CLFFT_COMPILER_NOT_AVAILABLE			= CL_COMPILER_NOT_AVAILABLE,
	CLFFT_DEVICE_NOT_AVAILABLE				= CL_DEVICE_NOT_AVAILABLE,
	CLFFT_DEVICE_NOT_FOUND					= CL_DEVICE_NOT_FOUND,
	CLFFT_SUCCESS							= CL_SUCCESS,
	//-------------------------- Extended status codes for clfft ----------------------------------------
	CLFFT_BUGCHECK =  4*1024,	/*!< Bugcheck. */
	CLFFT_NOTIMPLEMENTED,		/*!< Functionality is not implemented yet. */
	CLFFT_TRANSPOSED_NOTIMPLEMENTED, /*!< Transposed functionality is not implemented for this transformation. */
	CLFFT_FILE_NOT_FOUND,		/*!< Tried to open an existing file on the host system, but failed. */
	CLFFT_FILE_CREATE_FAILURE,	/*!< Tried to create a file on the host system, but failed. */
	CLFFT_VERSION_MISMATCH,		/*!< Version conflict between client and library. */
	CLFFT_INVALID_PLAN,			/*!< Requested plan could not be found. */
	CLFFT_DEVICE_NO_DOUBLE,		/*!< Double precision not supported on this device. */
	CLFFT_DEVICE_MISMATCH,		/*!< Attempt to run on a device using a plan baked for a different device. */
	CLFFT_ENDSTATUS				/* The last value of the enum, and marks the length of clfftStatus. */
};
typedef enum clfftStatus_ clfftStatus;

/*!  @brief The dimension of the input and output buffers that is fed into all FFT transforms */
typedef enum clfftDim_
{
	CLFFT_1D		= 1,		/*!< 1 Dimensional FFT transform (default). */
	CLFFT_2D,					/*!< 2 Dimensional FFT transform. */
	CLFFT_3D,					/*!< 3 Dimensional FFT transform. */
	ENDDIMENSION			/*!< The last value of the enum, and marks the length of clfftDim. */
} clfftDim;

/*!  @brief Specify the expected layouts of the buffers */
typedef enum clfftLayout_
{
	CLFFT_COMPLEX_INTERLEAVED	= 1,	/*!< An array of complex numbers, with real and imaginary components together (default). */
	CLFFT_COMPLEX_PLANAR,				/*!< Separate arrays of real components and imaginary components. */
	CLFFT_HERMITIAN_INTERLEAVED,		/*!< Compressed form of complex numbers; complex-conjugates are not stored, real and imaginary components are stored in the same array. */
	CLFFT_HERMITIAN_PLANAR,				/*!< Compressed form of complex numbers; complex-conjugates are not stored, real and imaginary components are stored in separate arrays. */
	CLFFT_REAL,							/*!< An array of real numbers, with no corresponding imaginary components. */
	ENDLAYOUT			/*!< The last value of the enum, and marks the length of clfftLayout. */
} clfftLayout;

/*!  @brief Specify the expected precision of each FFT.
 */
typedef enum clfftPrecision_
{
	CLFFT_SINGLE	= 1,	/*!< An array of complex numbers, with real and imaginary components saved as floats (default). */
	CLFFT_DOUBLE,			/*!< An array of complex numbers, with real and imaginary components saved as doubles. */
	CLFFT_SINGLE_FAST,		/*!< Faster implementation preferred. */
	CLFFT_DOUBLE_FAST,		/*!< Faster implementation preferred. */
	ENDPRECISION	/*!< The last value of the enum, and marks the length of clfftPrecision. */
} clfftPrecision;

/*!  @brief Specify the expected direction of each FFT, time or the frequency domains */
typedef enum clfftDirection_
{
	CLFFT_FORWARD	= -1,		/*!< FFT transform from time to frequency domain. */
	CLFFT_BACKWARD	= 1,		/*!< FFT transform from frequency to time domain. */
	CLFFT_MINUS		= -1,		/*!< Alias for the forward transform. */
	CLFFT_PLUS		= 1,		/*!< Alias for the backward transform. */
	ENDDIRECTION			/*!< The last value of the enum, and marks the length of clfftDirection. */
} clfftDirection;

/*!  @brief Specify wheter the input buffers are overwritten with results */
typedef enum clfftResultLocation_
{
	CLFFT_INPLACE		= 1,		/*!< Input and output buffers are the same (default). */
	CLFFT_OUTOFPLACE,				/*!< Input and output buffers are separate. */
	ENDPLACE				/*!< The last value of the enum, and marks the length of clfftPlaceness. */
} clfftResultLocation;

/*! @brief Determines whether the result is returned in original order. It is valid only for
dimensions greater than 1. */
typedef enum clfftResultTransposed_ {
	CLFFT_NOTRANSPOSE = 1,		/*!< The result is returned in the original order (default) */
	CLFFT_TRANSPOSED,			/*!< The result is transposed where transpose kernel is supported (possibly faster) */
	ENDTRANSPOSED			/*!< The last value of the enum, and marks the length of clfftResultTransposed */
} clfftResultTransposed;

/*! 	BitMasks to be used with clfftSetupData.debugFlags */
#define CLFFT_DUMP_PROGRAMS 0x1

/*! @brief Data structure that can be passed to clfftSetup() to control the behavior of the FFT runtime
 *  @details This structure contains values that can be initialized before instantiation of the FFT runtime
 *  with ::clfftSetup().  To initialize this structure, pass a pointer to a user struct to ::clfftInitSetupData( ),
 *  which clears the structure and sets the version member variables to the current values.
 */
struct clfftSetupData_
{
	cl_uint major;		/*!< Major version number of the project; signifies possible major API changes. */
	cl_uint minor;		/*!< Minor version number of the project; minor API changes that can break backward compatibility. */
	cl_uint patch;		/*!< Patch version number of the project; always incrementing number, signifies change over time. */

	/*! 	Bitwise flags that control the behavior of library debug logic. */
	cl_ulong debugFlags;  /*! This must be set to zero, except when debugging the clfft library.
	                       *  <p> debugFlags can be set to CLFFT_DUMP_PROGRAMS, in which case the dynamically generated OpenCL kernels are
	                       *  written to text files in the current working directory.  These files have a *.cl suffix.
	                       */
};
typedef struct clfftSetupData_ clfftSetupData;

/*! @brief Type of Callback function.
*/
typedef enum clfftCallbackType_
{
	PRECALLBACK,	/*!< Callback function is invoked only once for every point of input at the beginning of FFT transform. */
	POSTCALLBACK	/*!< Callback function is invoked only once for every point of output at the end of FFT transform. */
}clfftCallbackType;

/*!  @brief An abstract handle to the object that represents the state of the FFT(s) */
typedef size_t clfftPlanHandle;

#ifdef __cplusplus
extern "C" {
#endif
	/*! @brief Initialize a clfftSetupData struct for the client
	 *  @details clfftSetupData is passed to clfftSetup to control behavior of the FFT runtime.
	 *  @param[out] setupData Data structure is cleared and initialized with version information and default values
	 *  @return Enum describes the error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus clfftInitSetupData( clfftSetupData* setupData );


	/*! @brief Initialize the internal FFT resources.
	 *  @details The internal resources include FFT implementation caches kernels, programs, and buffers.
	 *  @param[in] setupData Data structure that is passed into the setup routine to control FFT generation behavior
	 * 	and debug functionality
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftSetup( const clfftSetupData* setupData );

	/*! @brief Release all internal resources.
	 *  @details Called when client is done with the FFT library, allowing the library to destroy all resources it has cached
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftTeardown( );

	/*! @brief Query the FFT library for version information
	 *  @details Returns the major, minor and patch version numbers associated with the FFT library
	 *  @param[out] major Major functionality change
	 *  @param[out] minor Minor functionality change
	 *  @param[out] patch Bug fixes, documentation changes, no new features introduced
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftGetVersion( cl_uint* major, cl_uint* minor, cl_uint* patch );

	/*! @brief Create a plan object initialized entirely with default values.
	 *  @details A plan is a repository of state for calculating FFT's.  Allows the runtime to pre-calculate kernels, programs
	 * 	and buffers and associate them with buffers of specified dimensions.
	 *  @param[out] plHandle Handle to the newly created plan
	 *  @param[in] context Client is responsible for providing an OpenCL context for the plan
	 *  @param[in] dim Dimensionality of the FFT transform; describes how many elements are in the array
	 *  @param[in] clLengths An array of length of size 'dim';  each array value describes the length of each dimension
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftCreateDefaultPlan( clfftPlanHandle* plHandle, cl_context context, const clfftDim dim,
								const size_t* clLengths );

	/*! @brief Create a copy of an existing plan.
	 *  @details This API allows a client to create a new plan based upon an existing plan.  This function can be used to
	 *  quickly create plans that are similar, but may differ slightly.
	 *  @param[out] out_plHandle Handle to the newly created plan that is based on in_plHandle
	 *  @param[in] new_context Client is responsible for providing a new context for the new plan
	 *  @param[in] in_plHandle Handle to a previously created plan that is to be copied
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftCopyPlan( clfftPlanHandle* out_plHandle, cl_context new_context, clfftPlanHandle in_plHandle );

	/*! @brief Prepare the plan for execution.
	 *  @details After all plan parameters are set, the client has the option of 'baking' the plan, which informs the runtime that
	 *  no more change to the parameters of the plan is expected, and the OpenCL kernels can be compiled.  This optional function
	 *  allows the client application to perform the OpenCL kernel compilation when the application is initialized instead of during the first
	 *  execution.
	 *  At this point, the clfft runtime applies all implimented optimizations, including
	 *  running kernel experiments on the devices in the plan context.
	 *  <p>  This function takes a long time to execute. If a plan is not baked before being executed,
	 *  the first call to clfftEnqueueTransform takes a long time to execute.
	 *  <p>  If any significant parameter of a plan is changed after the plan is baked (by a subsequent call to any one of
	 *  the functions that has the prefix "clfftSetPlan"), it is not considered an error.  Instead, the plan reverts back to
	 *  the unbaked state, discarding the benefits of the baking operation.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] numQueues Number of command queues in commQueueFFT; 0 is a valid value, in which case the client does not want
	 * 	the runtime to run load experiments and only pre-calculate state information
	 *  @param[in] commQueueFFT An array of cl_command_queues created by the client; the command queues must be a proper subset of
	 * 	the devices included in the plan context
	 *  @param[in] pfn_notify A function pointer to a notification routine. The notification routine is a callback function that
	 *  an application can register and is called when the program executable is built (successfully or unsuccessfully).
	 *  Currently, this parameter MUST be NULL or nullptr.
	 *  @param[in] user_data Passed as an argument when pfn_notify is called.
	 *  Currently, this parameter MUST be NULL or nullptr.
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftBakePlan( clfftPlanHandle plHandle, cl_uint numQueues, cl_command_queue* commQueueFFT,
							void (CL_CALLBACK *pfn_notify)(clfftPlanHandle plHandle, void *user_data), void* user_data );

	/*! @brief Release the resources of a plan.
	 *  @details A plan may include resources, such as kernels, programs, and buffers that consume memory.  When a plan
	 *  is no more needed, the client must release the plan.
	 *  @param[in,out] plHandle Handle to a previously created plan
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftDestroyPlan( clfftPlanHandle* plHandle );

	/*! @brief Retrieve the OpenCL context of a previously created plan.
	 *  @details The user must pass a reference to a cl_context variable, which is modified to point to a
	 *  context set in the specified plan.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] context Reference to the user allocated cl_context, which points to context set in the plan
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftGetPlanContext( const clfftPlanHandle plHandle, cl_context* context );

	/*! @brief Retrieve the floating point precision of the FFT data
	 *  @details The user must pass a reference to a clfftPrecision variable, which is set to the
	 *  precision of the FFT complex data in the plan.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] precision Reference to the user clfftPrecision enum
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftGetPlanPrecision( const clfftPlanHandle plHandle, clfftPrecision* precision );

	/*! @brief Set the floating point precision of the FFT data
	 *  @details Sets the floating point precision of the FFT complex data in the plan.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] precision Reference to the user clfftPrecision enum
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftSetPlanPrecision( clfftPlanHandle plHandle, clfftPrecision precision );

	/*! @brief Retrieve the scaling factor that is applied to the FFT data
	 *  @details The user must pass a reference to a cl_float variable, which is set to the
	 *  floating point scaling factor that is multiplied across the FFT data.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dir Direction of the applied scaling factor
	 *  @param[out] scale Reference to the user cl_float variable
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftGetPlanScale( const clfftPlanHandle plHandle, clfftDirection dir, cl_float* scale );

	/*! @brief Set the scaling factor that is applied to the FFT data
	 *  @details Sets the floating point scaling factor that is
	 *  multiplied across the FFT data.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dir Direction of the applied scaling factor
	 *  @param[in] scale Reference to the user cl_float variable
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftSetPlanScale( clfftPlanHandle plHandle, clfftDirection dir, cl_float scale );

	/*! @brief Retrieve the number of discrete arrays that the plan can concurrently handle
	 *  @details The user must pass a reference to a cl_uint variable, which is set to the
	 *  number of discrete arrays (1D or 2D) that is batched together for the plan
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] batchSize Number of discrete FFTs performed
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftGetPlanBatchSize( const clfftPlanHandle plHandle, size_t* batchSize );

	/*! @brief Set the number of discrete arrays that the plan can concurrently handle
	 *  @details Sets the plan property which sets the number of discrete arrays (1D or 2D)
	 *  that is batched together for the plan
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] batchSize Number of discrete FFTs performed
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftSetPlanBatchSize( clfftPlanHandle plHandle, size_t batchSize );

	/*! @brief Retrieve the dimensionality of the data that is transformed
	 *  @details Queries a plan object and retrieves the value of the dimensionality that the plan is set for.  A size is returned to
	 *  help the client allocate sufficient storage to hold the dimensions in a further call to clfftGetPlanLength
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] dim The dimensionality of the FFT to be transformed
	 *  @param[out] size Value to allocate an array to hold the FFT dimensions.
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftGetPlanDim( const clfftPlanHandle plHandle, clfftDim* dim, cl_uint* size );

	/*! @brief Set the dimensionality of the data that is transformed
	 *  @details Set the dimensionality of the data that is transformed by the plan
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dim The dimensionality of the FFT to be transformed
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftSetPlanDim( clfftPlanHandle plHandle, const clfftDim dim );

	/*! @brief Retrieve the length of each dimension of the FFT
	 *  @details The user must pass a reference to a size_t array, which is set to the
	 *  length of each discrete dimension of the FFT
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dim Dimension of the FFT; describes how many elements are in the clLengths array
	 *  @param[out] clLengths An array of length of size 'dim';  each array value describes the length of each dimension
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftGetPlanLength( const clfftPlanHandle plHandle, const clfftDim dim, size_t* clLengths );

	/*! @brief Set the length of each dimension of the FFT
	 *  @details Sets the plan property which is the length of each discrete dimension of the FFT
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dim The dimension of the FFT; describes how many elements are in the clLengths array
	 *  @param[in] clLengths An array of length of size 'dim';  each array value describes the length of each dimension
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftSetPlanLength( clfftPlanHandle plHandle, const clfftDim dim, const size_t* clLengths );

	/*! @brief Retrieve the distance between consecutive elements of input buffers in each dimension.
	 *  @details Depending on how the dimension is set in the plan (for 2D or 3D FFT), strideY or strideZ can be safely
	 *  ignored
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dim The dimension of the stride parameters; provides the number of elements in the array
	 *  @param[out] clStrides An array of strides, of size 'dim'.
	 */
	CLFFTAPI clfftStatus	clfftGetPlanInStride( const clfftPlanHandle plHandle, const clfftDim dim, size_t* clStrides );

	/*! @brief Set the distance between consecutive elements of input buffers in each dimension.
	 *  @details Set the plan properties which is the distance between elements in all dimensions of the input buffer
	 *  (units are in terms of clfftPrecision)
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dim The dimension of the stride parameters; provides the number of elements in the clStrides array
	 *  @param[in] clStrides An array of strides of size 'dim'. Usually, strideX=1 so that successive elements in the first dimension are stored contiguously.
	 * 	Typically, strideY=LenX and strideZ=LenX*LenY with the successive elements in the second and third dimensions stored in packed format.
	 *  See  @ref DistanceStridesandPitches for details.
	 */
	CLFFTAPI clfftStatus	clfftSetPlanInStride( clfftPlanHandle plHandle, const clfftDim dim, size_t* clStrides );

	/*! @brief Retrieve the distance between consecutive elements of output buffers in each dimension.
	 *  @details Depending on how the dimension is set in the plan (for 2D or 3D FFT), strideY or strideZ can be safely
	 *  ignored
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dim The dimension of the stride parameters; provides the number of elements in the clStrides array
	 *  @param[out] clStrides An array of strides, of size 'dim'.
	 */
	CLFFTAPI clfftStatus	clfftGetPlanOutStride( const clfftPlanHandle plHandle, const clfftDim dim, size_t* clStrides );

	/*! @brief Set the distance between consecutive elements of output buffers in a dimension.
	 *  @details Sets the plan properties which is the distance between elements in all dimensions of the output buffer
	 *  (units are in terms of clfftPrecision)
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dim The dimension of the stride parameters; provides the number of elements in the clStrides array
	 *  @param[in] clStrides An array of strides of size 'dim'.  Usually, strideX=1 so that successive elements in the first dimension are stored contiguously.
	 * 	Typically, strideY=LenX and strideZ=LenX*LenY cause the successive elements in the second and third dimensions be stored in packed format.
	 *  @sa clfftSetPlanInStride
	 */
	CLFFTAPI clfftStatus	clfftSetPlanOutStride( clfftPlanHandle plHandle, const clfftDim dim, size_t* clStrides );

	/*! @brief Retrieve the distance between array objects
	 *  @details Pitch is the distance between each discrete array object in an FFT array. This is only used
	 *  for 'array' dimensions in clfftDim; see clfftSetPlanDimension (units are in terms of clfftPrecision)
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] iDist The distance between the beginning elements of the discrete array objects in input buffer.
	 *  For contiguous arrays in memory, iDist=(strideX*strideY*strideZ)
	 *  @param[out] oDist The distance between the beginning elements of the discrete array objects in output buffer.
	 *  For contiguous arrays in memory, oDist=(strideX*strideY*strideZ)
	 */
	CLFFTAPI clfftStatus	clfftGetPlanDistance( const clfftPlanHandle plHandle, size_t* iDist, size_t* oDist );

	/*! @brief Set the distance between array objects
	 *  @details Pitch is the distance between each discrete array object in an FFT array. This is only used
	 *  for 'array' dimensions in clfftDim; see clfftSetPlanDimension (units are in terms of clfftPrecision)
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] iDist The distance between the beginning elements of the discrete array objects in input buffer.
	 *  For contiguous arrays in memory, iDist=(strideX*strideY*strideZ)
	 *  @param[out] oDist The distance between the beginning elements of the discrete array objects in output buffer.
	 *  For contiguous arrays in memory, oDist=(strideX*strideY*strideZ)
	 */
	CLFFTAPI clfftStatus	clfftSetPlanDistance( clfftPlanHandle plHandle, size_t iDist, size_t oDist );

	/*! @brief Retrieve the expected layout of the input and output buffers
	 *  @details Input and output buffers can be filled with either Hermitian, complex, or real numbers.  Complex numbers are stored
	 *  in various layouts; this function retrieves the layouts used by input and output
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] iLayout Indicates how the input buffers are laid out in memory
	 *  @param[out] oLayout Indicates how the output buffers are laid out in memory
	 */
	CLFFTAPI clfftStatus	clfftGetLayout( const clfftPlanHandle plHandle, clfftLayout* iLayout, clfftLayout* oLayout );

	/*! @brief Set the expected layout of the input and output buffers
	 *  @details Input and output buffers can be filled with either Hermitian, complex, or real numbers.  Complex numbers can be stored
	 *  in various layouts; this function informs the library what layouts to use for input and output
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] iLayout Indicates how the input buffers are laid out in memory
	 *  @param[in] oLayout Indicates how the output buffers are laid out in memory
	 */
	CLFFTAPI clfftStatus	clfftSetLayout( clfftPlanHandle plHandle, clfftLayout iLayout, clfftLayout oLayout );

	/*! @brief Retrieve whether the input buffers are to be overwritten with results
	 *  @details If the setting performs an in-place transform, the input buffers are overwritten with the results of the
	 *  transform.  If the setting performs an out-of-place transforms, the library looks for separate output buffers
	 *  on the Enqueue call.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] placeness Informs the library to either overwrite the input buffers with results or to write them in separate output buffers
	 */
	CLFFTAPI clfftStatus	clfftGetResultLocation( const clfftPlanHandle plHandle, clfftResultLocation* placeness );

	/*! @brief Set whether the input buffers are to be overwritten with results
	 *  @details If the setting performs an in-place transform, the input buffers are overwritten with the results of the
	 *  transform.  If the setting performs an out-of-place transforms, the library looks for separate output buffers
	 *  on the Enqueue call.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] placeness Informs the library to either overwrite the input buffers with results or to write them in separate output buffers
	 */
	CLFFTAPI clfftStatus	clfftSetResultLocation( clfftPlanHandle plHandle, clfftResultLocation placeness );

	/*! @brief Retrieve the final transpose setting of a multi-dimensional FFT
	 *  @details A multi-dimensional FFT transposes the data several times during calculation. If the client
	 *  does not care about the final transpose, to put data back in proper dimension, the final transpose can be skipped
	 *  to improve speed
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] transposed Specifies whether the final transpose can be skipped
	 */
	CLFFTAPI clfftStatus	clfftGetPlanTransposeResult( const clfftPlanHandle plHandle, clfftResultTransposed * transposed );

	/*! @brief Set the final transpose setting of a multi-dimensional FFT
	 *  @details A multi-dimensional FFT transposes the data several times during calculation.  If the client
	 *  does not care about the final transpose, to put data back in proper dimension, the final transpose can be skipped
	 *  to improve speed
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] transposed Specifies whether the final transpose can be skipped
	 */
	CLFFTAPI clfftStatus	clfftSetPlanTransposeResult( clfftPlanHandle plHandle, clfftResultTransposed transposed );


	/*! @brief Get buffer size (in bytes), which may be needed internally for an intermediate buffer
	 *  @details Very large FFT transforms may need multiple passes, and the operation needs a temporary buffer to hold
	 *  intermediate results. This function is only valid after the plan is baked, otherwise, an invalid operation error
	 *  is returned. If the returned buffersize is 0, the runtime needs no temporary buffer.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[out] buffersize Size in bytes for intermediate buffer
	 */
	CLFFTAPI clfftStatus clfftGetTmpBufSize( const clfftPlanHandle plHandle, size_t* buffersize );

	/*! @brief Register the callback parameters
	 *  @details Client can provide a callback function to do custom processing while reading input data and/or
	 *  writing output data. The callback function is provided as a string.
	 *  clFFT library incorporates the callback function string into the main FFT kernel. This function is used
	 *  by client to set the necessary parameters for callback
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] funcName Callback function name
	 *  @param[in] funcString Callback function in string form
	 *  @param[in] localMemSize Optional - Size (bytes) of the local memory used by callback function; pass 0 if no local memory is used
	 *  @param[in] callbackType Type of callback - Pre-Callback or Post-Callback
	 *  @param[in] userdata Supplementary data if any used by callback function
	 *  @param[in] numUserdataBuffers Number of userdata buffers
	 */
	CLFFTAPI clfftStatus clfftSetPlanCallback(clfftPlanHandle plHandle, const char* funcName, const char* funcString,
										int localMemSize, clfftCallbackType callbackType, cl_mem *userdata, int numUserdataBuffers);


	/*! @brief Enqueue an FFT transform operation, and return immediately (non-blocking)
	 *  @details This transform API function computes the FFT transform. It is non-blocking as it
	 *  only enqueues the OpenCL kernels for execution. The synchronization step must be managed by the user.
	 *  @param[in] plHandle Handle to a previously created plan
	 *  @param[in] dir Forward or backward transform
	 *  @param[in] numQueuesAndEvents Number of command queues in commQueues; number of expected events to be returned in outEvents
	 *  @param[in] commQueues An array of cl_command_queues created by the client; the command queues must be a proper subset of
	 * 	the devices included in the OpenCL context associated with the plan
	 *  @param[in] numWaitEvents Specify the number of elements in the eventWaitList array
	 *  @param[in] waitEvents Events for which the transform waits to complete before executing on the device
	 *  @param[out] outEvents The runtime fills this array with events corresponding one to one with the input command queues passed
	 *	in commQueues.  This parameter can have the value NULL or nullptr. When the value is NULL, the client is not interested in receiving notifications
	 *	when transforms are finished, otherwise, (if not NULL) the client is responsible for allocating this array with at least
	 *	as many elements as specified in numQueuesAndEvents.
	 *  @param[in] inputBuffers An array of cl_mem objects that contain data for processing by the FFT runtime. If the transform
	 *  is in-place, the FFT results overwrite the input buffers
	 *  @param[out] outputBuffers An array of cl_mem objects that store the results of out-of-place transforms. If the transform
	 *  is in-place, this parameter may be NULL or nullptr and is completely ignored
	 *  @param[in] tmpBuffer A cl_mem object that is reserved as a temporary buffer for FFT processing. If clTmpBuffers is NULL or nullptr,
	 *  and the library needs temporary storage, an internal temporary buffer is created on the fly managed by the library.
	 *  @return Enum describing error condition; superset of OpenCL error codes
	 */
	CLFFTAPI clfftStatus	clfftEnqueueTransform(
												clfftPlanHandle plHandle,
												clfftDirection dir,
												cl_uint numQueuesAndEvents,
												cl_command_queue* commQueues,
												cl_uint numWaitEvents,
												const cl_event* waitEvents,
												cl_event* outEvents,
												cl_mem* inputBuffers,
												cl_mem* outputBuffers,
												cl_mem tmpBuffer
												);

#ifdef __cplusplus
}
#endif

#endif
