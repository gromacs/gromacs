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


// clfft.lifetime.cpp : Functions that control the lifetime of the FFT library and their supporting functions
//

#include "stdafx.h"
#include "private.h"
#include "repo.h"
#include "../include/sharedLibrary.h"
#include "../statTimer/statisticalTimer.extern.h"

clfftStatus clfftInitSetupData( clfftSetupData* setupData )
{
	setupData->major	= clfftVersionMajor;
	setupData->minor	= clfftVersionMinor;
	setupData->patch	= clfftVersionPatch;
	setupData->debugFlags	= 0;

	return	CLFFT_SUCCESS;
}

//	Allow AMD's implementation of FFT's to allocate internal resources
clfftStatus	clfftSetup( const clfftSetupData* sData )
{
	//	Static data is not thread safe (to create), so we implement a lock to protect instantiation for the first call
	//	Implemented outside of FFTRepo::getInstance to minimize lock overhead; this is only necessary on first creation
	scopedLock sLock( FFTRepo::lockRepo(), _T( "FFTRepo::getInstance" ) );

	//	First invocation of this function will allocate the FFTRepo singleton; thereafter the object always exists
	FFTRepo& fftRepo	= FFTRepo::getInstance( );

	clfftInitRequestLibNoMemAlloc();
	clfftInitBinaryCache();

	//	Discover and load the timer module if present
	fftRepo.timerHandle = LoadSharedLibrary( "lib", "StatTimer", true );
	if( fftRepo.timerHandle )
	{
		//	Timer module discovered and loaded successfully
		//	Initialize function pointers to call into the shared module
		PFGETSTATTIMER pfGetStatTimer = reinterpret_cast< PFGETSTATTIMER > ( LoadFunctionAddr( fftRepo.timerHandle, "getStatTimer" ) );

		//	Create and initialize our timer class, if the external timer shared library loaded
		if( pfGetStatTimer )
		{
			fftRepo.pStatTimer = reinterpret_cast< GpuStatTimer* > ( pfGetStatTimer( CLFFT_GPU ) );
		}
	}

	// If the client has no setupData, we are done
	if( sData == NULL )
		return CLFFT_SUCCESS;

	//	Versioning checks commented out until necessary
	////	If the major version number between the client and library do not match, return mismatch
	//if( sData->major > clfftVersionMajor )
	//	return CLFFT_VERSION_MISMATCH;

	////	If the minor version number between the client and library do not match, return mismatch
	//if( sData->minor > clfftVersionMinor )
	//	return CLFFT_VERSION_MISMATCH;

	////	We ignore patch version number for version validation

	fftRepo.setupData	= *sData;

	return	CLFFT_SUCCESS;
}

//	Allow AMD's implementation of FFT's to destroy internal resources
clfftStatus	clfftTeardown( )
{
	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	fftRepo.releaseResources( );

	FreeSharedLibrary( fftRepo.timerHandle );

	return	CLFFT_SUCCESS;
}

clfftStatus clfftGetVersion( cl_uint* major, cl_uint* minor, cl_uint* patch )
{
	*major	= clfftVersionMajor;
	*minor	= clfftVersionMinor;
	*patch	= clfftVersionPatch;

	return	CLFFT_SUCCESS;
}
