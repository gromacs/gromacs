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


#pragma once
#if !defined( CLFFT_private_H )
#define CLFFT_private_H

#include <vector>
#include <string>
#include <locale>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <cassert>
#include "../include/clFFT.h"
#include "../include/unicode.compatibility.h"

#if defined(_MSC_VER)
	//	Microsoft Visual C++ compiler
	//
#define SPRINTF(_buffer, _count, _format,...)                        \
	_snprintf_s (_buffer, _count, _TRUNCATE, _format, __VA_ARGS__)
#elif defined( __GNUC__ )
	//	Gnu G++
	//
#define SPRINTF(_buffer, _count, _format,...)                   \
	{	size_t len = (_count)-1;                                \
		snprintf (_buffer, len, _format,__VA_ARGS__);           \
		_buffer[len] = 0;                                       \
	}
#else
#error Unknown/unsupported C++ compiler.
#endif

//	Creating a portable defintion of countof
//  This excludes mingw compilers; mingw32 does not have _countof
#if defined( _MSC_VER )
	#define countOf _countof
#else
	#define countOf( arr ) ( sizeof( arr ) / sizeof( arr[ 0 ] ) )
#endif

// This excludes mingw compilers; mingw32 does not have <intrin.h>
#if defined( _MSC_VER )
	#include <intrin.h>

	#if defined( _WIN64 )
		inline void BSF( unsigned long* index, size_t& mask )
		{
			_BitScanForward64( index, mask );
		}

		inline size_t AtomicAdd( volatile size_t* value, size_t op )
		{
			return _InterlockedExchangeAdd64( reinterpret_cast< volatile __int64* >( value ), op );
		}
	#else
		inline void BSF( unsigned long* index, size_t& mask )
		{
			_BitScanForward( index, mask );
		}

		inline size_t AtomicAdd( volatile size_t* value, size_t op )
		{
			return _InterlockedExchangeAdd( reinterpret_cast< volatile long* >( value ), op );
		}
	#endif
#elif defined( __GNUC__ )
	inline void BSF( unsigned long * index, size_t & mask )
	{
		*index = __builtin_ctz( mask );
	}

	inline size_t AtomicAdd( volatile size_t* value, size_t op )
	{
		return __sync_fetch_and_add( value, op );
	}
#endif

void clfftInitRequestLibNoMemAlloc();
bool clfftGetRequestLibNoMemAlloc();

void clfftInitBinaryCache();

//	This header file is not visible to clients, and contains internal structures and functions for use
//	by the FFT library.  Since this header is private to this implementation, there is no need to keep
//	strict C compliance.

//	Enum to help provide descriptive names to array indices, when indexing into our various vectors
enum clfftDim_Index
{
	DimX,				///< 1 Dimension
	DimY,				///< 2 Dimension
	DimZ,				///< 3 Dimension
	DimW,				///< 4th Dimension
	ENDDIMINDEX			///< This value will always be last, and marks the length of clfftDim_Index
};

template< typename FileStreamType, typename StringType >
class tofstreamRAII
{
	FileStreamType	outFile;
	StringType		fileName;

	public:
		tofstreamRAII( const StringType& name ): fileName( name )
		{
			outFile.open( fileName.c_str( ) );
		}

		~tofstreamRAII( )
		{
			outFile.close( );
		}

		StringType& getName( )
		{
			return fileName;
		}

		void setName( const StringType& name )
		{
			fileName = name;
		}

		FileStreamType& get( )
		{
			return outFile;
		}
};

//(currently) true if length is a power of 2,3,5,7,11,13
inline bool IsASupportedLength( size_t length )
{
	while( length > 1 )
	{
		if( length % 2 == 0 )
			length /= 2;
		else if( length % 3 == 0 )
			length /= 3;
		else if( length % 5 == 0 )
			length /= 5;
		else if( length % 7 == 0 )
			length /= 7;
		else if (length % 11 == 0)
			length /= 11;
		else if (length % 13 == 0)
			length /= 13;
		else
			return false;
	}
	return true;
}

inline tstring clfftErrorStatusAsString( const cl_int& status )
{
	switch( status )
	{
		case CLFFT_INVALID_GLOBAL_WORK_SIZE:
			return _T( "CLFFT_INVALID_GLOBAL_WORK_SIZE" );
		case CLFFT_INVALID_MIP_LEVEL:
			return _T( "CLFFT_INVALID_MIP_LEVEL" );
		case CLFFT_INVALID_BUFFER_SIZE:
			return _T( "CLFFT_INVALID_BUFFER_SIZE" );
		case CLFFT_INVALID_GL_OBJECT:
			return _T( "CLFFT_INVALID_GL_OBJECT" );
		case CLFFT_INVALID_OPERATION:
			return _T( "CLFFT_INVALID_OPERATION" );
		case CLFFT_INVALID_EVENT:
			return _T( "CLFFT_INVALID_EVENT" );
		case CLFFT_INVALID_EVENT_WAIT_LIST:
			return _T( "CLFFT_INVALID_EVENT_WAIT_LIST" );
		case CLFFT_INVALID_GLOBAL_OFFSET:
			return _T( "CLFFT_INVALID_GLOBAL_OFFSET" );
		case CLFFT_INVALID_WORK_ITEM_SIZE:
			return _T( "CLFFT_INVALID_WORK_ITEM_SIZE" );
		case CLFFT_INVALID_WORK_GROUP_SIZE:
			return _T( "CLFFT_INVALID_WORK_GROUP_SIZE" );
		case CLFFT_INVALID_WORK_DIMENSION:
			return _T( "CLFFT_INVALID_WORK_DIMENSION" );
		case CLFFT_INVALID_KERNEL_ARGS:
			return _T( "CLFFT_INVALID_KERNEL_ARGS" );
		case CLFFT_INVALID_ARG_SIZE:
			return _T( "CLFFT_INVALID_ARG_SIZE" );
		case CLFFT_INVALID_ARG_VALUE:
			return _T( "CLFFT_INVALID_ARG_VALUE" );
		case CLFFT_INVALID_ARG_INDEX:
			return _T( "CLFFT_INVALID_ARG_INDEX" );
		case CLFFT_INVALID_KERNEL:
			return _T( "CLFFT_INVALID_KERNEL" );
		case CLFFT_INVALID_KERNEL_DEFINITION:
			return _T( "CLFFT_INVALID_KERNEL_DEFINITION" );
		case CLFFT_INVALID_KERNEL_NAME:
			return _T( "CLFFT_INVALID_KERNEL_NAME" );
		case CLFFT_INVALID_PROGRAM_EXECUTABLE:
			return _T( "CLFFT_INVALID_PROGRAM_EXECUTABLE" );
		case CLFFT_INVALID_PROGRAM:
			return _T( "CLFFT_INVALID_PROGRAM" );
		case CLFFT_INVALID_BUILD_OPTIONS:
			return _T( "CLFFT_INVALID_BUILD_OPTIONS" );
		case CLFFT_INVALID_BINARY:
			return _T( "CLFFT_INVALID_BINARY" );
		case CLFFT_INVALID_SAMPLER:
			return _T( "CLFFT_INVALID_SAMPLER" );
		case CLFFT_INVALID_IMAGE_SIZE:
			return _T( "CLFFT_INVALID_IMAGE_SIZE" );
		case CLFFT_INVALID_IMAGE_FORMAT_DESCRIPTOR:
			return _T( "CLFFT_INVALID_IMAGE_FORMAT_DESCRIPTOR" );
		case CLFFT_INVALID_MEM_OBJECT:
			return _T( "CLFFT_INVALID_MEM_OBJECT" );
		case CLFFT_INVALID_HOST_PTR:
			return _T( "CLFFT_INVALID_HOST_PTR" );
		case CLFFT_INVALID_COMMAND_QUEUE:
			return _T( "CLFFT_INVALID_COMMAND_QUEUE" );
		case CLFFT_INVALID_QUEUE_PROPERTIES:
			return _T( "CLFFT_INVALID_QUEUE_PROPERTIES" );
		case CLFFT_INVALID_CONTEXT:
			return _T( "CLFFT_INVALID_CONTEXT" );
		case CLFFT_INVALID_DEVICE:
			return _T( "CLFFT_INVALID_DEVICE" );
		case CLFFT_INVALID_PLATFORM:
			return _T( "CLFFT_INVALID_PLATFORM" );
		case CLFFT_INVALID_DEVICE_TYPE:
			return _T( "CLFFT_INVALID_DEVICE_TYPE" );
		case CLFFT_INVALID_VALUE:
			return _T( "CLFFT_INVALID_VALUE" );
		case CLFFT_MAP_FAILURE:
			return _T( "CLFFT_MAP_FAILURE" );
		case CLFFT_BUILD_PROGRAM_FAILURE:
			return _T( "CLFFT_BUILD_PROGRAM_FAILURE" );
		case CLFFT_IMAGE_FORMAT_NOT_SUPPORTED:
			return _T( "CLFFT_IMAGE_FORMAT_NOT_SUPPORTED" );
		case CLFFT_IMAGE_FORMAT_MISMATCH:
			return _T( "CLFFT_IMAGE_FORMAT_MISMATCH" );
		case CLFFT_MEM_COPY_OVERLAP:
			return _T( "CLFFT_MEM_COPY_OVERLAP" );
		case CLFFT_PROFILING_INFO_NOT_AVAILABLE:
			return _T( "CLFFT_PROFILING_INFO_NOT_AVAILABLE" );
		case CLFFT_OUT_OF_HOST_MEMORY:
			return _T( "CLFFT_OUT_OF_HOST_MEMORY" );
		case CLFFT_OUT_OF_RESOURCES:
			return _T( "CLFFT_OUT_OF_RESOURCES" );
		case CLFFT_MEM_OBJECT_ALLOCATION_FAILURE:
			return _T( "CLFFT_MEM_OBJECT_ALLOCATION_FAILURE" );
		case CLFFT_COMPILER_NOT_AVAILABLE:
			return _T( "CLFFT_COMPILER_NOT_AVAILABLE" );
		case CLFFT_DEVICE_NOT_AVAILABLE:
			return _T( "CLFFT_DEVICE_NOT_AVAILABLE" );
		case CLFFT_DEVICE_NOT_FOUND:
			return _T( "CLFFT_DEVICE_NOT_FOUND" );
		case CLFFT_SUCCESS:
			return _T( "CLFFT_SUCCESS" );
		case CLFFT_NOTIMPLEMENTED:
			return _T( "CLFFT_NOTIMPLEMENTED" );
		case CLFFT_FILE_NOT_FOUND:
			return _T( "CLFFT_FILE_NOT_FOUND" );
		case CLFFT_FILE_CREATE_FAILURE:
			return _T( "CLFFT_FILE_CREATE_FAILURE" );
		case CLFFT_VERSION_MISMATCH:
			return _T( "CLFFT_VERSION_MISMATCH" );
		case CLFFT_INVALID_PLAN:
			return _T( "CLFFT_INVALID_PLAN" );
		default:
			return _T( "Error code not defined" );
		break;
	}
}

//	This is used to either wrap an OpenCL function call, or to explicitly check a variable for an OpenCL error condition.
//	If an error occurs, we issue a return statement to exit the calling function.
#if defined( _DEBUG )

#define OPENCL_V( fn, msg ) \
{ \
	clfftStatus vclStatus = static_cast< clfftStatus >( fn ); \
	switch( vclStatus ) \
	{ \
		case	CL_SUCCESS:		/**< No error */ \
			break; \
		default: \
		{ \
			terr << _T( "OPENCL_V< " ); \
			terr << clfftErrorStatusAsString( vclStatus ); \
			terr << _T( " > (" )<< static_cast<unsigned>( __LINE__ ) << _T( "): " ); \
			terr << msg << std::endl; \
			return	vclStatus; \
		} \
	} \
}

#else

#define OPENCL_V( fn, msg ) \
{ \
	clfftStatus vclStatus = static_cast< clfftStatus >( fn ); \
	switch( vclStatus ) \
	{ \
		case	CL_SUCCESS:		/**< No error */ \
			break; \
		default: \
		{ \
			return	vclStatus; \
		} \
	} \
}
#endif

static inline bool IsPo2 (size_t u) {
	return (u != 0) &&  (0 == (u & (u-1)));
}

template<typename T>
static inline T DivRoundingUp (T a, T b) {
	return (a + (b-1)) / b;
}

static inline size_t BitScanF (size_t n) {
	assert (n != 0);
	unsigned long tmp = 0;
	BSF (& tmp, n);
	return (size_t) tmp;
}

#define ARG_CHECK(_proposition)	\
{ bool btmp = (_proposition);	assert (btmp); if (! btmp)	return CLFFT_INVALID_ARG_VALUE; }

#define BUG_CHECK(_proposition)	\
	{ bool btmp = (_proposition);	assert (btmp); if (! btmp)	return CLFFT_BUGCHECK; }

#ifdef __cplusplus
extern "C" {
#endif

CLFFTAPI clfftStatus clfftLocalMemSize( const clfftPlanHandle plHandle, cl_ulong* local_mem_size );

/*! @brief Save to disk a file containing the contents of a baked plan.
*  @details A plan is a repository of state for calculating FFT's. Saves the details for a plan to allow the user
*	to easily recreate a plan and execute it without having to first build the kernel.
*  @param[in] plHandle Handle to the plan to be written to disk
*  @param[in] filename The desired name of the output file for the stored plan
*  @return Enum describing error condition; superset of OpenCL error codes
*/
CLFFTAPI clfftStatus	clfftWritePlanToDisk( clfftPlanHandle plHandle, const char* filename );

/*! @brief Read from disk a file containing the contents of a baked plan.
*  @details A plan is a repository of state for calculating FFT's. Reads the details for a plan from a file on disk and duplicates
*	the plan in the provided plan handle.
*  @param[out] plHandle Handle to the plan to be set to details from the file
*  @param[in] filename The name of the file containing the stored plan
*  @return Enum describing error condition; superset of OpenCL error codes
*/
CLFFTAPI clfftStatus	clfftReadPlanFromDisk( clfftPlanHandle plHandle, const char* filename );



#ifdef __cplusplus
}
#endif

#endif
