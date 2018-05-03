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


// clfft.repo.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "repo.h"

using std::map;
using std::string;

//	Static initialization of the plan count variable
size_t FFTRepo::planCount	= 1;

//	Handle/Address of the dynamic module that contains the timer, that we discover and load during runtime
void* FFTRepo::timerHandle	= NULL;
GpuStatTimer* FFTRepo::pStatTimer	= NULL;




clfftStatus FFTRepo::releaseResources( )
{
	scopedLock sLock( lockRepo(), _T( "releaseResources" ) );

	//	Release all handles to Kernels
	//
	for(Kernel_iterator iKern = mapKernels.begin( ); iKern != mapKernels.end( ); ++iKern )
	{
		cl_kernel k = iKern->second.kernel_fwd;
		iKern->second.kernel_fwd = NULL;
		if (NULL != k)
			clReleaseKernel( k );
		k = iKern->second.kernel_back;
		iKern->second.kernel_back = NULL;
		if (NULL != k)
			clReleaseKernel( k );

		if (NULL != iKern->second.kernel_fwd_lock)
		{
			delete iKern->second.kernel_fwd_lock;
			iKern->second.kernel_fwd_lock = NULL;
		}

		if (NULL != iKern->second.kernel_back_lock)
		{
			delete iKern->second.kernel_back_lock;
			iKern->second.kernel_back_lock = NULL;
		}
	}
	mapKernels.clear( );

	//	Release all handles to programs
	//
	for (fftRepo_iterator iProg = mapFFTs.begin( ); iProg != mapFFTs.end( ); ++iProg )
	{
        if (iProg->first.data != NULL)
        {
            const_cast<FFTRepoKey*>(&iProg->first)->deleteData();
        }

		cl_program p = iProg->second.clProgram;
		iProg->second.clProgram = NULL;
		if (NULL != p)
			clReleaseProgram (p);
	}

	//	Free all memory allocated in the repoPlans; represents cached plans that were not destroyed by the client
	//
	for( repoPlansType::iterator iter = repoPlans.begin( ); iter != repoPlans.end( ); ++iter )
	{
		FFTPlan* plan	= iter->second.first;
		lockRAII* lock	= iter->second.second;
		if( plan != NULL )
		{
			delete plan;
		}
		if( lock != NULL )
		{
			delete lock;
		}
	}

	//	Reset the plan count to zero because we are guaranteed to have destroyed all plans
	planCount	= 1;

	//	Release all strings
	mapFFTs.clear( );

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::setProgramCode( const clfftGenerators gen, const FFTKernelSignatureHeader * data, const std::string& kernel, const cl_device_id &device, const cl_context& planContext )
{
	scopedLock sLock( lockRepo(), _T( "setProgramCode" ) );

	FFTRepoKey key(gen, data, planContext, device);

    key.privatizeData();

	// Prefix copyright statement at the top of generated kernels
	std::stringstream ss;
	ss << 
		"/* ************************************************************************\n"
		" * Copyright 2013 Advanced Micro Devices, Inc.\n"
		" *\n"
		" * Licensed under the Apache License, Version 2.0 (the \"License\");\n"
		" * you may not use this file except in compliance with the License.\n"
		" * You may obtain a copy of the License at\n"
		" *\n"
		" * http://www.apache.org/licenses/LICENSE-2.0\n"
		" *\n"
		" * Unless required by applicable law or agreed to in writing, software\n"
		" * distributed under the License is distributed on an \"AS IS\" BASIS,\n"
		" * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.\n"
		" * See the License for the specific language governing permissions and\n"
		" * limitations under the License.\n"
		" * ************************************************************************/"
	<< std::endl << std::endl;

	std::string prefixCopyright = ss.str();

	fftRepoType::iterator it = mapFFTs.find(key);
	if (it == mapFFTs.end())
		mapFFTs[key].ProgramString = prefixCopyright + kernel;
	else
		key.deleteData();

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::getProgramCode( const clfftGenerators gen, const FFTKernelSignatureHeader * data, std::string& kernel, const cl_device_id &device, const cl_context& planContext )
{
	scopedLock sLock( lockRepo(), _T( "getProgramCode" ) );

	FFTRepoKey key(gen, data, planContext, device);

	fftRepo_iterator pos = mapFFTs.find( key);
	if( pos == mapFFTs.end( ) )
		return	CLFFT_FILE_NOT_FOUND;

  kernel = pos->second.ProgramString;
	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::setProgramEntryPoints( const clfftGenerators gen, const FFTKernelSignatureHeader * data,
	const char * kernel_fwd, const char * kernel_back, const cl_device_id &device, const cl_context& planContext  )
{
	scopedLock sLock( lockRepo(), _T( "setProgramEntryPoints" ) );

	FFTRepoKey key(gen, data, planContext, device);

	fftRepoValue& fft  = mapFFTs[ key ];
	fft.EntryPoint_fwd  = kernel_fwd;
	fft.EntryPoint_back = kernel_back;

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::getProgramEntryPoint( const clfftGenerators gen, const FFTKernelSignatureHeader * data,
			clfftDirection dir, std::string& kernel, const cl_device_id &device, const cl_context& planContext )
{
	scopedLock sLock( lockRepo(), _T( "getProgramEntryPoint" ) );

	FFTRepoKey key(gen, data, planContext, device);

	fftRepo_iterator pos = mapFFTs.find( key );
	if( pos == mapFFTs.end( ) )
		return	CLFFT_FILE_NOT_FOUND;

	switch (dir) {
	case CLFFT_FORWARD:
		kernel = pos->second.EntryPoint_fwd;
		break;
	case CLFFT_BACKWARD:
		kernel = pos->second.EntryPoint_back;
		break;
	default:
		assert (false);
		return CLFFT_INVALID_ARG_VALUE;
	}

	if (0 == kernel.size())
		return	CLFFT_FILE_NOT_FOUND;

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::setclProgram( const clfftGenerators gen, const FFTKernelSignatureHeader * data, const cl_program& prog, const cl_device_id &device, const cl_context& planContext )
{
	scopedLock sLock( lockRepo(), _T( "setclProgram" ) );

 	FFTRepoKey key(gen, data, planContext, device);

	fftRepo_iterator pos = mapFFTs.find( key );
	if( pos == mapFFTs.end( ) )
	{
		key.privatizeData(); // the key owns the data
		mapFFTs[ key ].clProgram = prog;
	}
	else {
		cl_program p = pos->second.clProgram;
		assert (NULL == p);
		if (NULL != p)
			clReleaseProgram (p);
		pos->second.clProgram = prog;
	}

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::getclProgram( const clfftGenerators gen, const FFTKernelSignatureHeader * data, cl_program& prog, const cl_device_id &device, const cl_context& planContext  )
{
	scopedLock sLock( lockRepo(), _T( "getclProgram" ) );

	FFTRepoKey key(gen, data, planContext, device);

	fftRepo_iterator pos = mapFFTs.find( key );
	if( pos == mapFFTs.end( ) )
		return	CLFFT_INVALID_PROGRAM;
	prog = pos->second.clProgram;
	if (NULL == prog)
		return	CLFFT_INVALID_PROGRAM;
  
  cl_context progContext;
  clGetProgramInfo(prog, CL_PROGRAM_CONTEXT, sizeof(cl_context), &progContext, NULL);
  if (planContext!=progContext)
    return	CLFFT_INVALID_PROGRAM;

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::setclKernel( cl_program prog, clfftDirection dir, const cl_kernel& kernel )
{
	scopedLock sLock( lockRepo(), _T( "setclKernel" ) );

	fftKernels & Kernels = mapKernels[ prog ];

	cl_kernel * pk;
	lockRAII ** kernelLock;

	switch (dir) {
	case CLFFT_FORWARD:
		pk = & Kernels.kernel_fwd;
		kernelLock = & Kernels.kernel_fwd_lock;
		break;
	case CLFFT_BACKWARD:
		pk = & Kernels.kernel_back;
		kernelLock = & Kernels.kernel_back_lock;
		break;
	default:
		assert (false);
		return CLFFT_INVALID_ARG_VALUE;
	}

	assert (NULL == *pk);
	if (NULL != *pk)
		clReleaseKernel( *pk );

	*pk = kernel;

	if (NULL != *kernelLock)
		 delete kernelLock;

	*kernelLock = new lockRAII;

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::getclKernel( cl_program prog, clfftDirection dir, cl_kernel& kernel, lockRAII*& kernelLock)
{
	scopedLock sLock( lockRepo(), _T( "getclKernel" ) );

	Kernel_iterator pos = mapKernels.find( prog );
	if (pos == mapKernels.end( ) )
		return	CLFFT_INVALID_KERNEL;

	switch (dir) {
	case CLFFT_FORWARD:
		kernel = pos->second.kernel_fwd;
		kernelLock = pos->second.kernel_fwd_lock;
		break;
	case CLFFT_BACKWARD:
		kernel = pos->second.kernel_back;
		kernelLock = pos->second.kernel_back_lock;
		break;
	default:
		assert (false);
		return CLFFT_INVALID_ARG_VALUE;
	}

	if (NULL == kernel)
		return	CLFFT_INVALID_KERNEL;

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::createPlan( clfftPlanHandle* plHandle, FFTPlan*& fftPlan )
{
	scopedLock sLock( lockRepo(), _T( "insertPlan" ) );

	//	We keep track of this memory in our own collection class, to make sure it's freed in releaseResources
	//	The lifetime of a plan is tracked by the client and is freed when the client calls ::clfftDestroyPlan()
	fftPlan	= new FFTPlan;

	//	We allocate a new lock here, and expect it to be freed in ::clfftDestroyPlan();
	//	The lifetime of the lock is the same as the lifetime of the plan
	lockRAII* lockPlan	= new lockRAII;

	//	Add and remember the fftPlan in our map
	repoPlans[ planCount ] = std::make_pair( fftPlan, lockPlan );

	//	Assign the user handle the plan count (unique identifier), and bump the count for the next plan
	*plHandle	= planCount++;

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::getPlan( clfftPlanHandle plHandle, FFTPlan*& fftPlan, lockRAII*& planLock )
{
	scopedLock sLock( lockRepo(), _T( "getPlan" ) );

	//	First, check if we have already created a plan with this exact same FFTPlan
	repoPlansType::iterator iter	= repoPlans.find( plHandle );
	if( iter == repoPlans.end( ) )
		return	CLFFT_INVALID_PLAN;

	//	If plan is valid, return fill out the output pointers
	fftPlan		= iter->second.first;
	planLock	= iter->second.second;

	return	CLFFT_SUCCESS;
}

clfftStatus FFTRepo::deletePlan( clfftPlanHandle* plHandle )
{
	scopedLock sLock( lockRepo(), _T( "deletePlan" ) );

	//	First, check if we have already created a plan with this exact same FFTPlan
	repoPlansType::iterator iter	= repoPlans.find( *plHandle );
	if( iter == repoPlans.end( ) )
		return	CLFFT_INVALID_PLAN;

	//	We lock the plan object while we are in the process of deleting it
	{
		scopedLock sLock( *iter->second.second, _T( "clfftDestroyPlan" ) );
		clReleaseContext( iter->second.first->context );

		//	Delete the FFTPlan
		delete iter->second.first;
	}

		//	Delete the lockRAII
	delete iter->second.second;

	//	Remove entry from our map object
	repoPlans.erase( iter );

	//	Clear the client's handle to signify that the plan is gone
	*plHandle = 0;

	return	CLFFT_SUCCESS;
}
