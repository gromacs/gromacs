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
#ifndef _STATISTICALTIMER_EXTERN_H_
#define _STATISTICALTIMER_EXTERN_H_
#include "../include/clFFT.h"
#include "statisticalTimer.h"

/**
 * \file clfft.StatisticalTimer.extern.h
 * \brief A timer class that provides a cross platform timer for use
 * in timing code progress with a high degree of accuracy.
 *	This class is implemented entirely in the header, to facilitate inclusion into multiple
 *	projects without needing to compile an object file for each project.
 */

// The following ifdef block is the standard way of creating macros which make exporting
// from a DLL simpler. All files within this DLL are compiled with the STATTIMER_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see
// STATTIMER_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#if defined( _WIN32 )
	#if !defined( __cplusplus )
		#define inline __inline
	#endif

    #if defined( CLFFT_STATIC )
        #define STATTIMER_API
    #elif defined( STATTIMER_EXPORTS )
        #define STATTIMER_API __declspec( dllexport )
    #else
        #define STATTIMER_API __declspec( dllimport )
    #endif
#else
	#define STATTIMER_API
#endif

//	The type of timer to be returned from ::getStatTimer( )
typedef enum clfftTimerType_
{
	CLFFT_GPU			= 1,
	CLFFT_CPU,
} clfftTimerType;

//	Table of typedef definitions for all exported functions from this shared module.
//	Clients of this module can use these typedefs to help create function pointers
//	that can be initialized to point to the functions exported from this module.
typedef baseStatTimer* (*PFGETSTATTIMER)( const clfftTimerType type );

	/**
	* \fn getInstance()
	* \brief This returns a reference to the singleton timer.  Guarantees only 1 timer class is ever
	*	instantiated within a compilable executable.
	*/
extern "C" STATTIMER_API baseStatTimer* getStatTimer( const clfftTimerType type );

#endif // _STATISTICALTIMER_EXTERN_H_
