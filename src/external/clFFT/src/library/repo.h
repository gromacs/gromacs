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
#if !defined( CLFFT_repo_H )
#define CLFFT_repo_H
#include <map>
#include "private.h"
#include "plan.h"
#include "lock.h"

#include "../statTimer/statisticalTimer.GPU.h"





//	This class contains objects that we wish to retain between individual calls into the FFT interface.
//	These objects will be shared across different individual FFT plans, and we wish to keep only one
//	copy of these programs, objects and events.  When the client decides that they either want to reset
//	the library or release all resources, this Repo will release all acquired resources and clean itself
//	up as much as it can.  It is implemented as a Singleton object.
class	FFTRepo
{

    struct FFTRepoKey
    {
        clfftGenerators gen;
        const FFTKernelSignatureHeader * data;
        cl_context context;
		cl_device_id device;
		bool dataIsPrivate;

        FFTRepoKey(clfftGenerators gen_, const FFTKernelSignatureHeader * data_, cl_context context_, cl_device_id device_)
            : gen(gen_), data(data_), context(context_), device(device_), dataIsPrivate(false)
        {

        }

        void privatizeData()
        {
            char * tmp = new char[data->datasize];
            ::memcpy(tmp, data, data->datasize);
            this->data = (FFTKernelSignatureHeader*) tmp;
			dataIsPrivate = true;
        }

        void deleteData()
        {
            if ( dataIsPrivate && (this->data != NULL) )
            {
				char *tmp = (char *)(this->data);
                delete[] tmp;
				this->data = 0;
            }            
        }

        bool operator<(const FFTRepoKey & b) const
        {
            const FFTRepoKey & a = *this;

            if (a.gen != b.gen)
            {
                return a.gen < b.gen;
            }
            else if (a.data->datasize != b.data->datasize)
            {
                return a.data->datasize < b.data->datasize;
            }
            else if (a.context != b.context)
            {
                return a.context < b.context;
            }
			else if (a.device != b.device)
			{
				return a.device < b.device;
			}
            else
            {
                return ::memcmp(a.data, b.data, a.data->datasize) < 0;
            }
        }
    };


	//	Structure containing all the data we need to remember for a specific invokation of a kernel
	//	generator
	struct fftRepoValue {
		std::string ProgramString;
		std::string EntryPoint_fwd;
		std::string EntryPoint_back;
		cl_program  clProgram;

		fftRepoValue ()
		:	clProgram (NULL)
		{}
	};

	//	Map structure to map parameters that a generator uses to a specific set of kernels that the generator
	//	has created
	//typedef std::pair< clfftGenerators, FFTKernelGenKeyParams > fftRepoKey;

	typedef std::map< FFTRepoKey, fftRepoValue > fftRepoType;
	typedef fftRepoType::iterator fftRepo_iterator;



	fftRepoType	mapFFTs;

	struct fftKernels {
		cl_kernel kernel_fwd;
		cl_kernel kernel_back;
		lockRAII* kernel_fwd_lock;
		lockRAII* kernel_back_lock;

		fftKernels ()
		:	kernel_fwd (NULL)
		,	kernel_back (NULL)
		,	kernel_fwd_lock(NULL)
		,	kernel_back_lock(NULL)
		{}
	};

	typedef std::map< cl_program, fftKernels > mapKernelType;
	typedef mapKernelType::iterator Kernel_iterator;
	mapKernelType mapKernels;

	//	All plans that the user creates over the course of using the library are stored here.
	//	Plans can be arbitrarily created and destroyed at anytime by the user, in arbitrary order, so vector
	//	does not seem appropriate, so a map was chosen because of the O(log N) search properties
	//	A lock object is created for each plan, such that any getter/setter can lock the 'plan' object before
	//	reading/writing its values.  The lock object is kept seperate from the plan object so that the lock
	//	object can be held the entire time a plan is getting destroyed in clfftDestroyPlan.
	typedef std::pair< FFTPlan*, lockRAII* > repoPlansValue;
	typedef std::map< clfftPlanHandle, repoPlansValue > repoPlansType;
	repoPlansType repoPlans;

	//	Static count of how many plans we have generated; always incrementing during the life of the library
	//	This is used as a unique identifier for plans
	static size_t planCount;

	// Private constructor to stop explicit instantiation
	FFTRepo( )
	{}

	// Private copy constructor to stop implicit instantiation
	FFTRepo( const FFTRepo& );

	// Private operator= to assure only 1 copy of singleton
	FFTRepo& operator=( const FFTRepo& );

	~FFTRepo( )
	{
		//	NOTE:  We can't release resources in our destructor because as a static object, the order of destruction of static objects
		//	is not guaranteed, and openCL might already have cleaned itself up.  When clFFT tries to free its resources, an access
		//	violation could occur.
		//releaseResources( );

		//	We should at least print out a warning message to the user if we are in our destructor and we still have resources
		//	bound.  This should give the user a clue to remember to call clfftTeardown( )
		if( (!mapKernels.empty( )) || (!mapFFTs.empty( )) )
		{
			terr << _T( "Warning:  Program terminating, but clFFT resources not freed." ) << std::endl;
			terr << _T( "Please consider explicitly calling clfftTeardown( )." ) << std::endl;
		}
	};

public:
	//	Used to make the FFTRepo struct thread safe; STL is not thread safe by default
	//	Optimally, we could use a lock object per STL struct, as two different STL structures
	//	can be modified at the same time, but a single lock object is easier and performance should
	//	still be good. This is implemented as a function returning a static local reference to
	//	assert that the lock must be instantiated before the result can be used.
	static lockRAII& lockRepo()
	{
		//	Static initialization of the repo lock variable
		static lockRAII lock(_T("FFTRepo"));
		return lock;
	}

	//	Our runtime library can instrument kernel timings with a GPU timer available in a shared module
	//	Handle/Address of the dynamic module that contains timers
	static void* timerHandle;

	//	Pointer to the timer class queried from the timer shared library
	static GpuStatTimer* pStatTimer;

	//	Global debug flags that the FFT engine can check to control various debug logic
	clfftSetupData setupData;

	//	Everybody who wants to access the Repo calls this function to get a repo reference
	static FFTRepo& getInstance( )
	{
		static FFTRepo fftRepo;
		return	fftRepo;
	};

	clfftStatus releaseResources( );

	clfftStatus setProgramCode( const clfftGenerators gen, const FFTKernelSignatureHeader * data, const std::string& kernel, const cl_device_id &device, const cl_context& planContext );
	clfftStatus getProgramCode( const clfftGenerators gen, const FFTKernelSignatureHeader * data, std::string& kernel, const cl_device_id &device, const cl_context& planContext );

	clfftStatus setProgramEntryPoints( const clfftGenerators gen, const FFTKernelSignatureHeader * data, const char * kernel_fwd, const char * kernel_back, const cl_device_id &device, const cl_context& planContext );
	clfftStatus getProgramEntryPoint( const clfftGenerators gen, const FFTKernelSignatureHeader * data, clfftDirection dir, std::string& kernel , const cl_device_id &device, const cl_context& planContext );

	clfftStatus setclProgram( const clfftGenerators gen, const FFTKernelSignatureHeader * data, const cl_program& prog, const cl_device_id &device, const cl_context& planContext );
	clfftStatus getclProgram( const clfftGenerators gen, const FFTKernelSignatureHeader * data, cl_program& prog, const cl_device_id &device, const cl_context& planContext );

	clfftStatus setclKernel ( cl_program prog, clfftDirection dir, const cl_kernel& kernel );
	clfftStatus getclKernel ( cl_program prog, clfftDirection dir, cl_kernel& kernel, lockRAII*& kernelLock);

	clfftStatus createPlan( clfftPlanHandle* plHandle, FFTPlan*& fftPlan );
	clfftStatus getPlan( clfftPlanHandle plHandle, FFTPlan*& fftPlan, lockRAII*& planLock );
	clfftStatus deletePlan( clfftPlanHandle* plHandle );
  

};

#endif

