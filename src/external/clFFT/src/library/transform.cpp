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


// clfft.transform.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "private.h"
#include "repo.h"
#include "plan.h"

//#define DEBUGGING

using std::vector;

clfftStatus clfftEnqueueTransform(
											clfftPlanHandle plHandle,
											clfftDirection dir,
											cl_uint numQueuesAndEvents,
											cl_command_queue* commQueues,
											cl_uint numWaitEvents,
											const cl_event* waitEvents,
											cl_event* outEvents,
											cl_mem* clInputBuffers,
											cl_mem* clOutputBuffers,
											cl_mem clTmpBuffers
											)
{
	cl_int status = CLFFT_SUCCESS;

	//	We do not currently support multiple command queues, which is necessary to support multi-gpu operations
	if( numQueuesAndEvents > 1 )
	{
		return CLFFT_NOTIMPLEMENTED;
	}

	FFTRepo& fftRepo	= FFTRepo::getInstance( );
	FFTPlan* fftPlan	= NULL;
	lockRAII* planLock	= NULL;

	//	At this point, the user wants to enqueue a plan to execute.  We lock the plan down now, such that
	//	after we finish baking the plan (if the user did not do that explicitely before), the plan cannot
	//	change again through the action of other thread before we enqueue this plan for execution.
	OPENCL_V( fftRepo.getPlan( plHandle, fftPlan, planLock ), _T( "fftRepo.getPlan failed" ) );
	scopedLock sLock( *planLock, _T( "clfftGetPlanBatchSize" ) );

	if( fftPlan->baked == false )
	{
		OPENCL_V( clfftBakePlan( plHandle, numQueuesAndEvents, commQueues, NULL, NULL ), _T( "Failed to bake plan" ) );
	}


	// get the device information
	cl_device_id q_device;
	clGetCommandQueueInfo(*commQueues, CL_QUEUE_DEVICE, sizeof(cl_device_id), &q_device, NULL);

	// verify if the current device is the same as the one used for baking the plan
	if(q_device != fftPlan->bakeDevice)
		return CLFFT_DEVICE_MISMATCH;


	if		(fftPlan->inputLayout == CLFFT_REAL)	dir = CLFFT_FORWARD;
	else if	(fftPlan->outputLayout == CLFFT_REAL)	dir = CLFFT_BACKWARD;


	// we do not check the user provided buffer at this release
	cl_mem localIntBuffer = clTmpBuffers;

	if( clTmpBuffers == NULL && fftPlan->tmpBufSize > 0 && fftPlan->intBuffer == NULL)
	{
		// create the intermediate buffers
		// The intermediate buffer is always interleave and packed
		// For outofplace operation, we have the choice not to create intermediate buffer
		// input ->(col+Transpose) output ->(col) output
		fftPlan->intBuffer = clCreateBuffer( fftPlan->context, CL_MEM_READ_WRITE,
			fftPlan->tmpBufSize, 0, &status);
		OPENCL_V( status, _T("Creating the intermediate buffer for large1D Failed") );
		fftPlan->libCreatedIntBuffer = true;

#if defined(DEBUGGING)
		std::cout << "One intermediate buffer is created" << std::endl;
#endif
	}

	if( localIntBuffer == NULL && fftPlan->intBuffer != NULL )
		localIntBuffer = fftPlan->intBuffer;

	if( fftPlan->intBufferRC == NULL && fftPlan->tmpBufSizeRC > 0 )
	{
		fftPlan->intBufferRC = clCreateBuffer( fftPlan->context, CL_MEM_READ_WRITE, fftPlan->tmpBufSizeRC, 0, &status);
		OPENCL_V( status, _T("Creating the intermediate buffer for large1D RC Failed") );
	}

	if( fftPlan->intBufferC2R == NULL && fftPlan->tmpBufSizeC2R > 0 )
	{
		fftPlan->intBufferC2R = clCreateBuffer( fftPlan->context, CL_MEM_READ_WRITE, fftPlan->tmpBufSizeC2R, 0, &status);
		OPENCL_V( status, _T("Creating the intermediate buffer for large1D YZ C2R Failed") );
	}

	//	The largest vector we can transform in a single pass
	//	depends on the GPU caps -- especially the amount of LDS
	//	available
	//
	size_t Large1DThreshold = 0;
	OPENCL_V(fftPlan->GetMax1DLength (&Large1DThreshold), _T("GetMax1DLength failed"));
	BUG_CHECK (Large1DThreshold > 1);

	//Large1DThreshold = 128;

	if(fftPlan->gen != Copy)
	switch( fftPlan->dim )
	{
		case CLFFT_1D:
		{
			if ( Is1DPossible(fftPlan->length[0], Large1DThreshold) )
				break;

			if( ( fftPlan->inputLayout == CLFFT_REAL ) && ( fftPlan->planTZ != 0) )
			{
					//First transpose
					// Input->tmp
					cl_event transTXOutEvents = NULL;
					OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
						waitEvents, &transTXOutEvents, clInputBuffers, &localIntBuffer, NULL ),
						_T("clfftEnqueueTransform for large1D transTX failed"));

					cl_mem *mybuffers;
					if (fftPlan->placeness==CLFFT_INPLACE)
						mybuffers = clInputBuffers;
					else
						mybuffers = clOutputBuffers;

#if defined(DEBUGGING)
								//  For debugging interleave data only,
								//  read the input buffer back into memory.
						clFinish(*commQueues);
								OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 0,
									NULL, NULL ),
									_T("Reading the result buffer failed") );
#endif

					//First Row
					//tmp->output
					cl_event rowXOutEvents = NULL;
					OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, 1,
						&transTXOutEvents, &rowXOutEvents, &localIntBuffer, &(fftPlan->intBufferRC), NULL ),
						_T("clfftEnqueueTransform for large1D rowX failed"));
					clReleaseEvent(transTXOutEvents);


#if defined(DEBUGGING)
								//  For debugging interleave data only,
								//  read the input buffer back into memory.
						clFinish(*commQueues);
								OPENCL_V( clEnqueueReadBuffer( *commQueues, *mybuffers, CL_TRUE, 0, 536870912, &temp[ 0 ], 0,
									NULL, NULL ),
									_T("Reading the result buffer failed") );
#endif

					//Second Transpose
					// output->tmp
					cl_event transTYOutEvents = NULL;
					OPENCL_V( clfftEnqueueTransform( fftPlan->planTY, dir, numQueuesAndEvents, commQueues, 1,
						&rowXOutEvents, &transTYOutEvents, &(fftPlan->intBufferRC), &localIntBuffer, NULL ),
						_T("clfftEnqueueTransform for large1D transTY failed"));
					clReleaseEvent(rowXOutEvents);


#if defined(DEBUGGING)
								//  For debugging interleave data only,
								//  read the input buffer back into memory.
						clFinish(*commQueues);
								OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 0,
									NULL, NULL ),
									_T("Reading the result buffer failed") );
#endif

					//Second Row
					//tmp->tmp, inplace
					cl_event rowYOutEvents = NULL;
					OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1,
						&transTYOutEvents, &rowYOutEvents, &localIntBuffer, &(fftPlan->intBufferRC), NULL ),
						_T("clfftEnqueueTransform for large1D rowY failed"));
					clReleaseEvent(transTYOutEvents);

#if defined(DEBUGGING)
								//  For debugging interleave data only,
								//  read the input buffer back into memory.
						clFinish(*commQueues);
								OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 0,
									NULL, NULL ),
									_T("Reading the result buffer failed") );
#endif

					//Third Transpose
					// tmp->output
					OPENCL_V( clfftEnqueueTransform( fftPlan->planTZ, dir, numQueuesAndEvents, commQueues, 1,
						&rowYOutEvents, outEvents, &(fftPlan->intBufferRC), mybuffers, NULL ),
						_T("clfftEnqueueTransform for large1D transTZ failed"));
					clReleaseEvent(rowYOutEvents);
			}
			else if ( fftPlan->inputLayout == CLFFT_REAL )
			{
				cl_event colOutEvents = NULL;
				cl_event copyInEvents = NULL;

				// First pass
				// column with twiddle first, OUTOFPLACE, + transpose
				OPENCL_V( clfftEnqueueTransform( fftPlan->planX, CLFFT_FORWARD, numQueuesAndEvents, commQueues, numWaitEvents,
					waitEvents, &colOutEvents, clInputBuffers, &(fftPlan->intBufferRC), localIntBuffer),
					_T("clfftEnqueueTransform large1D col pass failed"));


				cl_mem *out_local;
				out_local = (fftPlan->placeness==CLFFT_INPLACE) ? clInputBuffers : clOutputBuffers;


				// another column FFT output, INPLACE
				OPENCL_V(clfftEnqueueTransform(fftPlan->planY, CLFFT_FORWARD, numQueuesAndEvents, commQueues, 1, &colOutEvents,
					&copyInEvents, &(fftPlan->intBufferRC), &(fftPlan->intBufferRC), localIntBuffer),
					_T("clfftEnqueueTransform large1D second column failed"));
				clReleaseEvent(colOutEvents);

				// copy from full complex to hermitian
				OPENCL_V(clfftEnqueueTransform(fftPlan->planRCcopy, CLFFT_FORWARD, numQueuesAndEvents, commQueues, 1, &copyInEvents,
					outEvents, &(fftPlan->intBufferRC), out_local, localIntBuffer),
					_T("clfftEnqueueTransform large1D RC copy failed"));
				clReleaseEvent(copyInEvents);


			}
			else if( fftPlan->outputLayout == CLFFT_REAL )
			{
				cl_event colOutEvents = NULL;
				cl_event copyOutEvents = NULL;

				if (fftPlan->planRCcopy)
				{
					// copy from hermitian to full complex
					OPENCL_V(clfftEnqueueTransform(fftPlan->planRCcopy, CLFFT_BACKWARD, numQueuesAndEvents, commQueues, numWaitEvents,
						waitEvents, &copyOutEvents, clInputBuffers, &(fftPlan->intBufferRC), localIntBuffer),
						_T("clfftEnqueueTransform large1D RC copy failed"));

					// First pass
					// column with twiddle first, INPLACE,
					OPENCL_V(clfftEnqueueTransform(fftPlan->planX, CLFFT_BACKWARD, numQueuesAndEvents, commQueues, 1,
						&copyOutEvents, &colOutEvents, &(fftPlan->intBufferRC), &(fftPlan->intBufferRC), localIntBuffer),
						_T("clfftEnqueueTransform large1D col pass failed"));
					clReleaseEvent(copyOutEvents);
				}
				else
				{
					// First pass
					// column with twiddle first, INPLACE,
					OPENCL_V(clfftEnqueueTransform(fftPlan->planX, CLFFT_BACKWARD, numQueuesAndEvents, commQueues, numWaitEvents,
						waitEvents, &colOutEvents, clInputBuffers, &(fftPlan->intBufferRC), localIntBuffer),
						_T("clfftEnqueueTransform large1D col pass failed"));
					clReleaseEvent(copyOutEvents);
				}

				cl_mem *out_local;
				out_local = (fftPlan->placeness==CLFFT_INPLACE) ? clInputBuffers : clOutputBuffers;

				// another column FFT output, OUTOFPLACE + transpose
				OPENCL_V( clfftEnqueueTransform( fftPlan->planY, CLFFT_BACKWARD, numQueuesAndEvents, commQueues, 1, &colOutEvents,
					outEvents, &(fftPlan->intBufferRC), out_local, localIntBuffer ),
					_T("clfftEnqueueTransform large1D second column failed"));
				clReleaseEvent(colOutEvents);

			}
			else
			{
#if defined(DEBUGGING)
				// For debugging interleave data only, initialize the intermediate buffer
				// to a data pattern.  This will show which data in the buffer
				// are being written by the kernel
				//
				size_t buffSizeBytes_complex = fftPlan->tmpBufSize;
				size_t buffersize = buffSizeBytes_complex/sizeof( std::complex< float > );
				std::vector<std::complex< float> > temp(buffersize);

				for (size_t u = 0; u < buffersize; ++u) {
					temp[u] = std::complex<float> (float(u+1), float(buffersize-u));
				}

				if (fftPlan->large1D == 0)
				{
					//First time usage, we can initialize tmp buffer
					OPENCL_V(clEnqueueWriteBuffer( *commQueues,
						localIntBuffer,
						CL_TRUE,		// blocking write
						0,
						buffSizeBytes_complex,
						&temp[0],
						0,
						NULL,
						NULL), _T("clEnqueueWriteBuffer failed") );
				}
#endif

				if (fftPlan->transflag)
				{
					//First transpose
					// Input->tmp
					cl_event transTXOutEvents = NULL;
					if(fftPlan->allOpsInplace)
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
							waitEvents, &transTXOutEvents, clInputBuffers, NULL, NULL ),
							_T("clfftEnqueueTransform for large1D transTX failed"));
					}
					else
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
							waitEvents, &transTXOutEvents, clInputBuffers, &localIntBuffer, NULL ),
							_T("clfftEnqueueTransform for large1D transTX failed"));
					}

					cl_mem *mybuffers;
					if (fftPlan->placeness==CLFFT_INPLACE)
						mybuffers = clInputBuffers;
					else
						mybuffers = clOutputBuffers;

#if defined(DEBUGGING)
								//  For debugging interleave data only,
								//  read the input buffer back into memory.
						clFinish(*commQueues);
								OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 0,
									NULL, NULL ),
									_T("Reading the result buffer failed") );
#endif

					//First Row
					//tmp->output
					cl_event rowXOutEvents = NULL;
					if(fftPlan->allOpsInplace)
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, 1,
							&transTXOutEvents, &rowXOutEvents, clInputBuffers, NULL, NULL ),
							_T("clfftEnqueueTransform for large1D rowX failed"));
					}
					else
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, 1,
							&transTXOutEvents, &rowXOutEvents, &localIntBuffer, mybuffers, NULL ),
							_T("clfftEnqueueTransform for large1D rowX failed"));
					}
					clReleaseEvent(transTXOutEvents);


#if defined(DEBUGGING)
								//  For debugging interleave data only,
								//  read the input buffer back into memory.
						clFinish(*commQueues);
								OPENCL_V( clEnqueueReadBuffer( *commQueues, *mybuffers, CL_TRUE, 0, 536870912, &temp[ 0 ], 0,
									NULL, NULL ),
									_T("Reading the result buffer failed") );
#endif

					//Second Transpose
					// output->tmp
					cl_event transTYOutEvents = NULL;
					if(fftPlan->allOpsInplace)
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planTY, dir, numQueuesAndEvents, commQueues, 1,
							&rowXOutEvents, &transTYOutEvents, clInputBuffers, NULL, NULL ),
							_T("clfftEnqueueTransform for large1D transTY failed"));
					}
					else
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planTY, dir, numQueuesAndEvents, commQueues, 1,
							&rowXOutEvents, &transTYOutEvents, mybuffers, &localIntBuffer, NULL ),
							_T("clfftEnqueueTransform for large1D transTY failed"));
					}
					clReleaseEvent(rowXOutEvents);


#if defined(DEBUGGING)
								//  For debugging interleave data only,
								//  read the input buffer back into memory.
						clFinish(*commQueues);
								OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 0,
									NULL, NULL ),
									_T("Reading the result buffer failed") );
#endif

					//Second Row
					//tmp->tmp, inplace
					cl_event rowYOutEvents = NULL;
					if(fftPlan->allOpsInplace)
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1,
							&transTYOutEvents, &rowYOutEvents, clInputBuffers, NULL, NULL ),
							_T("clfftEnqueueTransform for large1D rowY failed"));
					}
					else
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1,
							&transTYOutEvents, &rowYOutEvents, &localIntBuffer, NULL, NULL ),
							_T("clfftEnqueueTransform for large1D rowY failed"));
					}
					clReleaseEvent(transTYOutEvents);

#if defined(DEBUGGING)
								//  For debugging interleave data only,
								//  read the input buffer back into memory.
						clFinish(*commQueues);
								OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 0,
									NULL, NULL ),
									_T("Reading the result buffer failed") );
#endif

					//Third Transpose
					// tmp->output
					if(fftPlan->allOpsInplace)
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planTZ, dir, numQueuesAndEvents, commQueues, 1,
							&rowYOutEvents, outEvents, clInputBuffers, NULL, NULL ),
							_T("clfftEnqueueTransform for large1D transTZ failed"));
					}
					else
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planTZ, dir, numQueuesAndEvents, commQueues, 1,
							&rowYOutEvents, outEvents, &localIntBuffer, mybuffers, NULL ),
							_T("clfftEnqueueTransform for large1D transTZ failed"));
					}
					clReleaseEvent(rowYOutEvents);

				}
				else
				{
					if (fftPlan->large1D == 0)
					{
						if(fftPlan->planCopy)
						{
							// Transpose OUTOFPLACE
							cl_event transTXOutEvents = NULL;
							OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
								waitEvents, &transTXOutEvents, clInputBuffers, &localIntBuffer, NULL ),
								_T("clfftEnqueueTransform for large1D transTX failed"));

#if defined(DEBUGGING)
									//  For debugging interleave data only,
									//  read the input buffer back into memory.
							clFinish(*commQueues);
									OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 0,
										NULL, NULL ),
										_T("Reading the result buffer failed") );
#endif

							// FFT INPLACE
							cl_event rowXOutEvents = NULL;
							OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, 1,
								&transTXOutEvents, &rowXOutEvents, &localIntBuffer, NULL, NULL),
								_T("clfftEnqueueTransform large1D first row pass failed"));
							clReleaseEvent(transTXOutEvents);

#if defined(DEBUGGING)
									//  For debugging interleave data only,
									//  read the input buffer back into memory.
							clFinish(*commQueues);
									OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 0,
										NULL, NULL ),
										_T("Reading the result buffer failed") );
#endif

							// FFT INPLACE
							cl_event colYOutEvents = NULL;
							OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &rowXOutEvents,
								&colYOutEvents, &localIntBuffer, NULL, NULL ),
								_T("clfftEnqueueTransform large1D second column failed"));
							clReleaseEvent(rowXOutEvents);
									
#if defined(DEBUGGING)
									//  For debugging interleave data only,
									//  read the input buffer back into memory.
							clFinish(*commQueues);
									OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 0,
										NULL, NULL ),
										_T("Reading the result buffer failed") );
#endif

							cl_mem *mybuffers;
							if (fftPlan->placeness==CLFFT_INPLACE)
								mybuffers = clInputBuffers;
							else
								mybuffers = clOutputBuffers;
						
							// Copy kernel
							OPENCL_V( clfftEnqueueTransform( fftPlan->planCopy, dir, numQueuesAndEvents, commQueues, 1, &colYOutEvents,
								outEvents, &localIntBuffer, mybuffers, NULL ),
								_T("clfftEnqueueTransform large1D copy failed"));
							clReleaseEvent(colYOutEvents);
						}
						else
						{
							cl_event colOutEvents = NULL;
							// First pass
							// column with twiddle first, OUTOFPLACE
							OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
								waitEvents, &colOutEvents, clInputBuffers, &localIntBuffer, NULL),
								_T("clfftEnqueueTransform large1D col pass failed"));

#if defined(DEBUGGING)
							// debug purpose, interleave input <-> interleave output
							// read the intermediate buffer and print part of it.
							OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 1,
								&colOutEvents, NULL ),
								_T("Reading the result buffer failed") );
#endif
							if(fftPlan->planTZ)
							{
								cl_event rowYOutEvents = NULL;
								OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &colOutEvents,
									&rowYOutEvents, &localIntBuffer, NULL, NULL ),
									_T("clfftEnqueueTransform large1D second row failed"));

								if (fftPlan->placeness == CLFFT_INPLACE)
								{
									OPENCL_V( clfftEnqueueTransform( fftPlan->planTZ, dir, numQueuesAndEvents, commQueues, 1, &rowYOutEvents,
										outEvents, &localIntBuffer, clInputBuffers, NULL ),
										_T("clfftEnqueueTransform large1D trans3 failed"));
								}
								else
								{
									OPENCL_V( clfftEnqueueTransform( fftPlan->planTZ, dir, numQueuesAndEvents, commQueues, 1, &rowYOutEvents,
										outEvents, &localIntBuffer, clOutputBuffers, NULL ),
										_T("clfftEnqueueTransform large1D trans3 failed"));
								}
						
								clReleaseEvent(rowYOutEvents);

							}
							else
							{
								//another column FFT output, OUTOFPLACE + transpose
								if (fftPlan->placeness == CLFFT_INPLACE)
								{
									OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &colOutEvents,
										outEvents, &localIntBuffer, clInputBuffers, NULL ),
										_T("clfftEnqueueTransform large1D second column failed"));

#if defined(DEBUGGING)
									//  For debugging interleave data only,
									//  read the input buffer back into memory.
									OPENCL_V( clEnqueueReadBuffer( *commQueues, clInputBuffers[0], CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 1,
										outEvents, NULL ),
										_T("Reading the result buffer failed") );
#endif
								}
								else
								{
#if defined(DEBUGGING)
								// debug purpose, interleave input <-> interleave output
								OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 1,
									&colOutEvents, NULL ),
									_T("Reading the result buffer failed") );
#endif
									OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &colOutEvents,
										outEvents, &localIntBuffer, clOutputBuffers, NULL ),
										_T("clfftEnqueueTransform large1D second column failed"));

#if defined(DEBUGGING)
									//  For debugging interleave data only, read back the output buffer
									//
									OPENCL_V( clEnqueueReadBuffer( *commQueues, clOutputBuffers[0], CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 1,
										outEvents, NULL ),
										_T("Reading the result buffer failed") );
#endif
								}
							}

							clReleaseEvent(colOutEvents);
						}
					}
					else
					{
						cl_event colOutEvents = NULL;

						// second pass for huge 1D
						// column with twiddle first, OUTOFPLACE, + transpose
						OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
							waitEvents, &colOutEvents, &localIntBuffer, clOutputBuffers, localIntBuffer),
							_T("clfftEnqueueTransform Huge1D col pass failed"));
#if defined(DEBUGGING)
						// debug purpose, interleave input <-> interleave output
						OPENCL_V( clEnqueueReadBuffer( *commQueues, clOutputBuffers[0], CL_TRUE, 0, buffSizeBytes_complex, &temp[ 0 ], 1,
							&colOutEvents, NULL ),
							_T("Reading the result buffer failed") );
#endif

						OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &colOutEvents,
							outEvents, clOutputBuffers, clOutputBuffers, localIntBuffer ),
							_T("clfftEnqueueTransform large1D second column failed"));

						clReleaseEvent(colOutEvents);
					}
				}
			}

			if( fftRepo.pStatTimer )
			{
				fftRepo.pStatTimer->AddSample( plHandle, fftPlan, NULL, 0, NULL, std::vector< size_t >( ), std::vector< size_t >() );
			}

			return	CLFFT_SUCCESS;

		}
		case CLFFT_2D:
		{
			// if transpose kernel, we will fall below
			if (fftPlan->transflag && !(fftPlan->planTX)) break;

			if ( (fftPlan->gen == Transpose_NONSQUARE ) &&
				 (fftPlan->nonSquareKernelType == NON_SQUARE_TRANS_PARENT) )
			{
				cl_event stage1OutEvents = NULL;

				OPENCL_V(clfftEnqueueTransform(fftPlan->planTX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
					waitEvents, &stage1OutEvents, clInputBuffers, NULL, NULL),
					_T("clfftEnqueueTransform stage1 failed"));

				OPENCL_V(clfftEnqueueTransform(fftPlan->planTY, dir, numQueuesAndEvents, commQueues, 1,
					&stage1OutEvents, outEvents, clInputBuffers, NULL, NULL),
					_T("clfftEnqueueTransform stage2 failed"));
				clReleaseEvent(stage1OutEvents);

				if (fftRepo.pStatTimer)
				{
					fftRepo.pStatTimer->AddSample(plHandle, fftPlan, NULL, 0, NULL, std::vector< size_t >(), std::vector< size_t >());
				}

				return	CLFFT_SUCCESS;
			}

			cl_event rowOutEvents = NULL;

#if defined(DEBUGGING)
			size_t buffersize = fftPlan->length[0] * fftPlan->length[1] * fftPlan->batchsize;
			if (fftPlan->length.size() > 2) buffersize *= fftPlan->length[2];
			//size_t buffSizeBytes=sizeof( std::complex< float > )*buffersize;
			//std::vector< std::complex< float > > output2( buffersize );
			size_t buffSizeBytes=sizeof( float) * buffersize;
			//std::vector<float> output2(buffersize*2);
			float *output2 = new float[buffersize*2];
#endif
#if defined(DEBUGGING)
			OPENCL_V( clEnqueueReadBuffer( *commQueues, clInputBuffers[0], CL_TRUE, 0, buffSizeBytes, &output2[ 0 ], 0,
				NULL, NULL ),
				_T("Reading the result buffer failed") );

			if (fftPlan->placeness == CLFFT_OUTOFPLACE)
			{
				OPENCL_V( clEnqueueReadBuffer( *commQueues, clOutputBuffers[0], CL_TRUE, 0, buffSizeBytes, &output2[ 0 ], 0,
					NULL, NULL ),
					_T("Reading the result buffer failed") );
			}
#endif
			if (fftPlan->transflag)
			{//first time set up transpose kernel for 2D
				//First row
				OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
					waitEvents, &rowOutEvents, clInputBuffers, clOutputBuffers, NULL ),
					_T("clfftEnqueueTransform for row failed"));

				cl_mem *mybuffers;

				if (fftPlan->placeness==CLFFT_INPLACE)
					mybuffers = clInputBuffers;
				else
					mybuffers = clOutputBuffers;

#if defined(DEBUGGING)
				OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[ 0 ], 0,
					NULL, NULL ),
					_T("Reading the result buffer failed") );
#endif

				cl_event transXOutEvents = NULL;
				cl_event colOutEvents = NULL;

				if (!fftPlan->transpose_in_2d_inplace)
				{
					//First transpose
					OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
						&transXOutEvents, mybuffers, &localIntBuffer, NULL ),
						_T("clfftEnqueueTransform for first transpose failed"));
					clReleaseEvent(rowOutEvents);

#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
#endif

					if (fftPlan->transposed == CLFFT_NOTRANSPOSE)
					{
						//Second Row transform
						OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
							&colOutEvents, &localIntBuffer, NULL, NULL ),
							_T("clfftEnqueueTransform for second row failed"));
						clReleaseEvent(transXOutEvents);

#if defined(DEBUGGING)
						OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
							NULL, NULL ),
							_T("Reading the result buffer failed") );
#endif

						//Second transpose
						OPENCL_V( clfftEnqueueTransform( fftPlan->planTY, dir, numQueuesAndEvents, commQueues, 1, &colOutEvents,
							outEvents, &localIntBuffer, mybuffers, NULL ),
							_T("clfftEnqueueTransform for second transpose failed"));
						clReleaseEvent(colOutEvents);

#if defined(DEBUGGING)
						OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
							NULL, NULL ),
							_T("Reading the result buffer failed") );
#endif
					}
					else
					{
						//Second Row transform
						OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
							outEvents, &localIntBuffer, mybuffers, NULL ),
							_T("clfftEnqueueTransform for second row failed"));
						clReleaseEvent(transXOutEvents);
					}
				}
				else
				{
					// First Transpose
					OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
						&transXOutEvents, mybuffers, NULL, NULL ),
						_T("clfftEnqueueTransform for first transpose failed"));
					clReleaseEvent(rowOutEvents);

					if (fftPlan->transposed == CLFFT_NOTRANSPOSE)
					{
						//Second Row transform
						OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
							&colOutEvents, mybuffers, NULL, NULL ),
							_T("clfftEnqueueTransform for Second Row failed"));
						clReleaseEvent(transXOutEvents);

						//Second transpose
						OPENCL_V( clfftEnqueueTransform( fftPlan->planTY, dir, numQueuesAndEvents, commQueues, 1, &colOutEvents,
							outEvents, mybuffers, NULL, NULL ),
							_T("clfftEnqueueTransform for second transpose failed"));
						clReleaseEvent(colOutEvents);
					}
					else
					{
						//Second Row transform
						OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
							outEvents, mybuffers, NULL, NULL ),
							_T("clfftEnqueueTransform for second row failed"));
						clReleaseEvent(transXOutEvents);
					}

				}
			}
			else
			{

				if ( (fftPlan->large2D || fftPlan->length.size()>2) &&
					(fftPlan->inputLayout != CLFFT_REAL) && (fftPlan->outputLayout != CLFFT_REAL))
				{
					if (fftPlan->placeness==CLFFT_INPLACE)
					{
						//deal with row first
						OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
							waitEvents, &rowOutEvents, clInputBuffers, NULL, localIntBuffer ),
							_T("clfftEnqueueTransform for row failed"));

						//deal with column
						OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
							outEvents, clInputBuffers, NULL, localIntBuffer ),
							_T("clfftEnqueueTransform for column failed"));
					}
					else
					{
						//deal with row first
						OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
							waitEvents, &rowOutEvents, clInputBuffers, clOutputBuffers, localIntBuffer ),
							_T("clfftEnqueueTransform for row failed"));

						//deal with column
						OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
							outEvents, clOutputBuffers, NULL, localIntBuffer ),
							_T("clfftEnqueueTransform for column failed"));

					}
				}
				else
				{
					if(fftPlan->inputLayout == CLFFT_REAL)
					{
						if(fftPlan->planTX)
						{
							//First row
							OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
								waitEvents, &rowOutEvents, clInputBuffers, clOutputBuffers, NULL ),
								_T("clfftEnqueueTransform for row failed"));

							cl_mem *mybuffers;

							if (fftPlan->placeness==CLFFT_INPLACE)
								mybuffers = clInputBuffers;
							else
								mybuffers = clOutputBuffers;

#if defined(DEBUGGING)
							OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[ 0 ], 0,
								NULL, NULL ),
								_T("Reading the result buffer failed") );
#endif

							cl_event transXOutEvents = NULL;
							cl_event colOutEvents = NULL;


							//First transpose
							OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
								&transXOutEvents, mybuffers, &localIntBuffer, NULL ),
								_T("clfftEnqueueTransform for first transpose failed"));
							// clReleaseEvent(rowOutEvents);

#if defined(DEBUGGING)
							OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
								NULL, NULL ),
								_T("Reading the result buffer failed") );
#endif


							//Second Row transform
							OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
								&colOutEvents, &localIntBuffer, NULL, NULL ),
								_T("clfftEnqueueTransform for second row failed"));
							clReleaseEvent(transXOutEvents);

#if defined(DEBUGGING)
							OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
								NULL, NULL ),
								_T("Reading the result buffer failed") );
#endif

							//Second transpose
							OPENCL_V( clfftEnqueueTransform( fftPlan->planTY, dir, numQueuesAndEvents, commQueues, 1, &colOutEvents,
								outEvents, &localIntBuffer, mybuffers, NULL ),
								_T("clfftEnqueueTransform for second transpose failed"));
							clReleaseEvent(colOutEvents);

#if defined(DEBUGGING)
							OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
								NULL, NULL ),
								_T("Reading the result buffer failed") );
#endif

						}
						else
						{
							if (fftPlan->placeness==CLFFT_INPLACE)
							{
								// deal with row
								OPENCL_V( clfftEnqueueTransform( fftPlan->planX, CLFFT_FORWARD, numQueuesAndEvents, commQueues, numWaitEvents,
									waitEvents, &rowOutEvents, clInputBuffers, NULL, localIntBuffer ),
									_T("clfftEnqueueTransform for row failed"));

								// deal with column
								OPENCL_V( clfftEnqueueTransform( fftPlan->planY, CLFFT_FORWARD, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
									outEvents, clInputBuffers, NULL, localIntBuffer ),
									_T("clfftEnqueueTransform for column failed"));
							}
							else
							{
								// deal with row
								OPENCL_V( clfftEnqueueTransform( fftPlan->planX, CLFFT_FORWARD, numQueuesAndEvents, commQueues, numWaitEvents,
									waitEvents, &rowOutEvents, clInputBuffers, clOutputBuffers, localIntBuffer ),
									_T("clfftEnqueueTransform for row failed"));

								// deal with column
								OPENCL_V( clfftEnqueueTransform( fftPlan->planY, CLFFT_FORWARD, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
									outEvents, clOutputBuffers, NULL, localIntBuffer ),
									_T("clfftEnqueueTransform for column failed"));
							}
						}
					}
					else if(fftPlan->outputLayout == CLFFT_REAL)
					{
						if(fftPlan->planTY)
						{
							cl_mem *mybuffers;

							if ( (fftPlan->placeness==CLFFT_INPLACE) ||
								 ((fftPlan->placeness==CLFFT_OUTOFPLACE) && (fftPlan->length.size() > 2)) )
								mybuffers = clInputBuffers;
							else
								mybuffers = &(fftPlan->intBufferC2R);

							cl_event transYOutEvents = NULL;
							cl_event transXOutEvents = NULL;

							//First transpose
							OPENCL_V( clfftEnqueueTransform( fftPlan->planTY, dir, numQueuesAndEvents, commQueues, numWaitEvents, 
								waitEvents, &transYOutEvents, clInputBuffers, &localIntBuffer, NULL ),
								_T("clfftEnqueueTransform for first transpose failed"));
					

#if defined(DEBUGGING)
							OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
								NULL, NULL ),
								_T("Reading the result buffer failed") );
#endif

							//First row
							OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &transYOutEvents, 
								&rowOutEvents, &localIntBuffer, NULL, NULL ),
								_T("clfftEnqueueTransform for col failed"));
							clReleaseEvent(transYOutEvents);


#if defined(DEBUGGING)
							OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[ 0 ], 0,
								NULL, NULL ),
								_T("Reading the result buffer failed") );
#endif

							//Second transpose
							OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
								&transXOutEvents, &localIntBuffer, mybuffers, NULL ),
								_T("clfftEnqueueTransform for second transpose failed"));
							

#if defined(DEBUGGING)
							OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
								NULL, NULL ),
								_T("Reading the result buffer failed") );
#endif


							//Second Row transform
							if(fftPlan->placeness == CLFFT_INPLACE)
							{
								OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
									outEvents, clInputBuffers, NULL, NULL ),
									_T("clfftEnqueueTransform for second row failed"));
							}
							else
							{
								OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
									outEvents, mybuffers, clOutputBuffers, NULL ),
									_T("clfftEnqueueTransform for second row failed"));
							}
							clReleaseEvent(transXOutEvents);
#if defined(DEBUGGING)
							OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
								NULL, NULL ),
								_T("Reading the result buffer failed") );
#endif


						}
						else
						{
							cl_mem *out_local, *int_local, *out_y;

							if(fftPlan->placeness == CLFFT_INPLACE)
							{
								out_local = NULL;
								int_local = NULL;
								out_y = clInputBuffers;
							}
							else
							{
								if(fftPlan->length.size() > 2)
								{
									out_local = clOutputBuffers;
									int_local = NULL;
									out_y = clInputBuffers;
								}
								else
								{
									out_local = clOutputBuffers;
									int_local = &(fftPlan->intBufferC2R);
									out_y = int_local;
								}
							}


							// deal with column
							OPENCL_V( clfftEnqueueTransform( fftPlan->planY, CLFFT_BACKWARD, numQueuesAndEvents, commQueues, numWaitEvents,
								waitEvents, &rowOutEvents, clInputBuffers, int_local, localIntBuffer ),
								_T("clfftEnqueueTransform for row failed"));

							// deal with row
							OPENCL_V( clfftEnqueueTransform( fftPlan->planX, CLFFT_BACKWARD, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
								outEvents, out_y, out_local, localIntBuffer ),
								_T("clfftEnqueueTransform for column failed"));
						}

					}
					else
					{
						//deal with row first
						OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
							waitEvents, &rowOutEvents, clInputBuffers, &localIntBuffer, NULL ),
							_T("clfftEnqueueTransform for row failed"));


						if (fftPlan->placeness==CLFFT_INPLACE)
						{
							//deal with column
							OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
								outEvents, &localIntBuffer, clInputBuffers, NULL ),
								_T("clfftEnqueueTransform for column failed"));
						}
						else
						{
							//deal with column
							OPENCL_V( clfftEnqueueTransform( fftPlan->planY, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
								outEvents, &localIntBuffer, clOutputBuffers, NULL ),
								_T("clfftEnqueueTransform for column failed"));

			#if defined(DEBUGGING)
							OPENCL_V( clEnqueueReadBuffer( *commQueues, clOutputBuffers[0], CL_TRUE, 0, buffSizeBytes, &output2[ 0 ], 1,
								outEvents, NULL ),
								_T("Reading the result buffer failed") );
			#endif
						}
					}
				}

				clReleaseEvent(rowOutEvents);

			}


			if( fftRepo.pStatTimer )
			{
				fftRepo.pStatTimer->AddSample( plHandle, fftPlan, NULL, 0, NULL, std::vector< size_t >( ), std::vector< size_t >() );
			}

			return	CLFFT_SUCCESS;
		}
		case CLFFT_3D:
		{
			cl_event rowOutEvents = NULL;

#if defined(DEBUGGING)
			size_t buffersize = fftPlan->length[0] * fftPlan->length[1] *fftPlan->length[2] *fftPlan->batchsize;
			size_t buffSizeBytes=sizeof( std::complex< float > )*buffersize;
			std::vector< std::complex< float > > output3( buffersize );
#endif
			if(fftPlan->inputLayout == CLFFT_REAL)
			{
				if(fftPlan->planTX)
				{
					//First row
					OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
						waitEvents, &rowOutEvents, clInputBuffers, clOutputBuffers, localIntBuffer ),
						_T("clfftEnqueueTransform for row failed"));

					cl_mem *mybuffers;

					if (fftPlan->placeness==CLFFT_INPLACE)
						mybuffers = clInputBuffers;
					else
						mybuffers = clOutputBuffers;

#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[ 0 ], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
#endif

					cl_event transXOutEvents = NULL;
					cl_event colOutEvents = NULL;


					//First transpose
					OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
						&transXOutEvents, mybuffers, &localIntBuffer, NULL ),
						_T("clfftEnqueueTransform for first transpose failed"));
					// clReleaseEvent(rowOutEvents);

#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
#endif


					//Second Row transform
					OPENCL_V( clfftEnqueueTransform( fftPlan->planZ, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
						&colOutEvents, &localIntBuffer, NULL, NULL ),
						_T("clfftEnqueueTransform for second row failed"));
					clReleaseEvent(transXOutEvents);

#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
#endif

					//Second transpose
					OPENCL_V( clfftEnqueueTransform( fftPlan->planTY, dir, numQueuesAndEvents, commQueues, 1, &colOutEvents,
						outEvents, &localIntBuffer, mybuffers, NULL ),
						_T("clfftEnqueueTransform for second transpose failed"));
					clReleaseEvent(colOutEvents);

#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
#endif

				}
				else
				{
					cl_mem *tmp_local, *out_local;

					tmp_local = (fftPlan->placeness==CLFFT_INPLACE) ? NULL : clOutputBuffers;
					out_local = (fftPlan->placeness==CLFFT_INPLACE) ? clInputBuffers : clOutputBuffers;

					//deal with 2D row first
					OPENCL_V( clfftEnqueueTransform( fftPlan->planX, CLFFT_FORWARD, numQueuesAndEvents, commQueues, numWaitEvents,
						waitEvents, &rowOutEvents, clInputBuffers, tmp_local, localIntBuffer ),
						_T("clfftEnqueueTransform for 3D-XY row failed"));

					//deal with 1D Z column
					OPENCL_V( clfftEnqueueTransform( fftPlan->planZ, CLFFT_FORWARD, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
						outEvents, out_local, NULL, localIntBuffer ),
						_T("clfftEnqueueTransform for 3D-Z column failed"));
				}

			}
			else if(fftPlan->outputLayout == CLFFT_REAL)
			{
				if(fftPlan->planTZ)
				{
					cl_mem *mybuffers;

					if (fftPlan->placeness==CLFFT_INPLACE)
						mybuffers = clInputBuffers;
					else
						mybuffers = &(fftPlan->intBufferC2R);

					cl_event transZOutEvents = NULL;
					cl_event transXOutEvents = NULL;

					//First transpose
					OPENCL_V( clfftEnqueueTransform( fftPlan->planTZ, dir, numQueuesAndEvents, commQueues, numWaitEvents, 
						waitEvents, &transZOutEvents, clInputBuffers, &localIntBuffer, NULL ),
						_T("clfftEnqueueTransform for first transpose failed"));
					

#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
#endif

					//First row
					OPENCL_V( clfftEnqueueTransform( fftPlan->planZ, dir, numQueuesAndEvents, commQueues, 1, &transZOutEvents, 
						&rowOutEvents, &localIntBuffer, NULL, NULL ),
						_T("clfftEnqueueTransform for col failed"));
					clReleaseEvent(transZOutEvents);


#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, mybuffers[0], CL_TRUE, 0, buffSizeBytes*2, &output2[ 0 ], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
#endif

					//Second transpose
					OPENCL_V( clfftEnqueueTransform( fftPlan->planTX, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
						&transXOutEvents, &localIntBuffer, mybuffers, NULL ),
						_T("clfftEnqueueTransform for second transpose failed"));
							

#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
#endif


					//Second Row transform
					if(fftPlan->placeness == CLFFT_INPLACE)
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
							outEvents, clInputBuffers, NULL, NULL ),
							_T("clfftEnqueueTransform for second row failed"));
					}
					else
					{
						OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, 1, &transXOutEvents,
							outEvents, mybuffers, clOutputBuffers, NULL ),
							_T("clfftEnqueueTransform for second row failed"));
					}
					clReleaseEvent(transXOutEvents);
#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, localIntBuffer, CL_TRUE, 0, buffSizeBytes*2, &output2[0], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
#endif


				}
				else
				{
					cl_mem *out_local, *int_local, *out_z;

					if(fftPlan->placeness == CLFFT_INPLACE)
					{
						out_local = NULL;
						int_local = NULL;
						out_z = clInputBuffers;
					}
					else
					{
						out_local = clOutputBuffers;
						int_local = &(fftPlan->intBufferC2R);
						out_z = int_local;
					}

					//deal with 1D Z column first
					OPENCL_V( clfftEnqueueTransform( fftPlan->planZ, CLFFT_BACKWARD, numQueuesAndEvents, commQueues, numWaitEvents,
						waitEvents, &rowOutEvents, clInputBuffers, int_local, localIntBuffer ),
						_T("clfftEnqueueTransform for 3D-Z column failed"));

					//deal with 2D row
					OPENCL_V( clfftEnqueueTransform( fftPlan->planX, CLFFT_BACKWARD, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
						outEvents, out_z, out_local, localIntBuffer ),
						_T("clfftEnqueueTransform for 3D-XY row failed"));
				}
			}
			else
			{
				if (fftPlan->placeness==CLFFT_INPLACE)
				{
					//deal with 2D row first
					OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
						waitEvents, &rowOutEvents, clInputBuffers, NULL, localIntBuffer ),
						_T("clfftEnqueueTransform for 3D-XY row failed"));

					//deal with 1D Z column
					OPENCL_V( clfftEnqueueTransform( fftPlan->planZ, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
						outEvents, clInputBuffers, NULL, localIntBuffer ),
						_T("clfftEnqueueTransform for 3D-Z column failed"));
				}
				else
				{
	#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, clOutputBuffers[0], CL_TRUE, 0, buffSizeBytes, &output3[ 0 ], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
	#endif

					//deal with 2D row first
					OPENCL_V( clfftEnqueueTransform( fftPlan->planX, dir, numQueuesAndEvents, commQueues, numWaitEvents,
						waitEvents, &rowOutEvents, clInputBuffers, clOutputBuffers, localIntBuffer ),
						_T("clfftEnqueueTransform for 3D-XY row failed"));

	#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, clOutputBuffers[0], CL_TRUE, 0, buffSizeBytes, &output3[ 0 ], 0,
						NULL, NULL ),
						_T("Reading the result buffer failed") );
	#endif

					//deal with 1D Z column
					OPENCL_V( clfftEnqueueTransform( fftPlan->planZ, dir, numQueuesAndEvents, commQueues, 1, &rowOutEvents,
						outEvents, clOutputBuffers, NULL, localIntBuffer ),
						_T("clfftEnqueueTransform for 3D-Z column failed"));
	#if defined(DEBUGGING)
					OPENCL_V( clEnqueueReadBuffer( *commQueues, clOutputBuffers[0], CL_TRUE, 0, buffSizeBytes, &output3[ 0 ], 1,
						outEvents, NULL ),
						_T("Reading the result buffer failed") );
	#endif
				}
			}

			clReleaseEvent(rowOutEvents);

			if( fftRepo.pStatTimer )
			{
				fftRepo.pStatTimer->AddSample( plHandle, fftPlan, NULL, 0, NULL, std::vector< size_t >( ), std::vector< size_t >() );
			}

			return	CLFFT_SUCCESS;
		}
	}

	return fftPlan->action->enqueue(plHandle,
                                        dir,
                                        numQueuesAndEvents,
                                        commQueues,
                                        numWaitEvents,
                                        waitEvents,
                                        outEvents,
                                        clInputBuffers,
                                        clOutputBuffers);
}
