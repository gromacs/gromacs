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
#ifndef _STATISTICALTIMER_GPU_H_
#define _STATISTICALTIMER_GPU_H_
#include <iosfwd>
#include <vector>
#include <algorithm>
#include <cmath>
#include "statisticalTimer.h"
#include "../library/plan.h"

/**
 * \file clfft.StatisticalTimer.GPU.h
 * \brief A timer class that provides a cross platform timer for use
 * in timing code progress with a high degree of accuracy.
 *	This class is implemented entirely in the header, to facilitate inclusion into multiple
 *	projects without needing to compile an object file for each project.
 */

struct StatData
{
	cl_kernel kernel;
	cl_ulong deltaNanoSec;
	double doubleNanoSec;
	size_t batchSize;
	clfftDim dim;
	clfftPlanHandle plHandle;
	clfftPlanHandle planX;
	clfftPlanHandle planY;
	clfftPlanHandle planZ;
	clfftPlanHandle planTX;
	clfftPlanHandle planTY;
	clfftPlanHandle planTZ;

	clfftPlanHandle planRCcopy;
	clfftPlanHandle planCopy;

	clfftGenerators gen;

	std::vector< size_t > lengths;
	std::vector< size_t > inStride;
	std::vector< size_t > outStride;
	size_t iDist;
	size_t oDist;
	clfftResultLocation placeness;
	std::vector< size_t > enqueueLocalWorkSize;
	std::vector< size_t > enqueueWorkSize;
	std::vector< cl_event > outEvents;

	StatData( ): deltaNanoSec( 0 )
	{}

	StatData( clfftPlanHandle id, FFTPlan* plan, cl_kernel kern, cl_uint nEv, cl_event* Ev,
		const std::vector< size_t >& gWorkSize, const std::vector< size_t >& lWorkSize):
		deltaNanoSec( 0 ), kernel( kern ), batchSize( plan->batchsize ), dim( plan->dim ),
		plHandle( id ), planX( plan->planX ), planY( plan->planY ), planZ( plan->planZ ),
		planTX( plan->planTX ), planTY( plan->planTY ), planTZ( plan->planTZ ),
		planRCcopy( plan->planRCcopy ), planCopy( plan->planCopy ), gen(plan->gen),
		inStride( plan->inStride ), outStride( plan->outStride ), iDist( plan->iDist ), oDist( plan->oDist ),
		lengths( plan->length ), enqueueWorkSize( gWorkSize ), enqueueLocalWorkSize( lWorkSize ), placeness( plan->placeness )
	{
		for( cl_uint e = 0; e < nEv; ++e )
		{
			outEvents.push_back( Ev[ e ] );
		}
	}

	double calcFlops( )
	{
		size_t	fftLength = 0;
		size_t	dimIndex = 0;

		if( dim == CLFFT_1D )
		{
			fftLength	= lengths.at( 0 );
			dimIndex	= 1;
		}
		else if( dim == CLFFT_2D )
		{
			fftLength	= lengths.at( 0 ) * lengths.at( 1 );
			dimIndex	= 2;
		}
		else if( dim == CLFFT_3D )
		{
			fftLength	= lengths.at( 0 ) * lengths.at( 1 ) * lengths.at( 2 );
			dimIndex	= 3;
		}

		size_t cumulativeBatch = 1;
		for( ; dimIndex < lengths.size(); ++dimIndex )
		{
			cumulativeBatch *= std::max< size_t >( 1, lengths[ dimIndex ] );
		}
		cumulativeBatch *= batchSize;

		double flops	= cumulativeBatch * 5 * fftLength * ( log( static_cast< double >( fftLength ) ) / log( 2.0 ) );

		return flops;
	}

};

//	Sorting operator for struct StatData, such that it can be used in a map
bool operator<( const StatData& lhs, const StatData& rhs);

class GpuStatTimer : public baseStatTimer
{
	//	Typedefs to handle the data that we store
	typedef std::vector< StatData > StatDataVec;
	typedef std::vector< StatDataVec > PerEnqueueVec;

	//	In order to calculate statistics <std. dev.>, we need to keep a history of our timings
	std::vector< PerEnqueueVec > timerData;

	//	Typedefs to handle the identifiers we use for our timers
	typedef	std::pair< std::string, cl_uint > idPair;
	typedef	std::vector< idPair > idVector;
	idVector labelID;

	//	Between each Start/Stop pair, we need to count how many AddSamples were made.
	size_t currSample, currRecord;

	//	Saved sizes for our vectors, used in Reset() to reallocate vectors
	StatDataVec::size_type	nEvents, nSamples;
	size_t currID;

	/**
	 * \fn GpuStatTimer()
	 * \brief Constructor for StatisticalTimer that initializes the class
	 *	This is private so that user code cannot create their own instantiation.  Instead, you
	 *	must go through getInstance( ) to get a reference to the class.
	 */
	GpuStatTimer( );

	/**
	 * \fn ~GpuStatTimer()
	 * \brief Destructor for StatisticalTimer that cleans up the class
	 */
	~GpuStatTimer( );

	/**
	 * \fn GpuStatTimer(const StatisticalTimer& )
	 * \brief Copy constructors do not make sense for a singleton, disallow copies
	 */
	GpuStatTimer( const GpuStatTimer& );

	/**
	 * \fn operator=( const StatisticalTimer& )
	 * \brief Assignment operator does not make sense for a singleton, disallow assignments
	 */
	GpuStatTimer& operator=( const GpuStatTimer& );

	friend std::ostream& operator<<( std::ostream& os, const GpuStatTimer& s );

	//	Calculate the average/mean of data for a given event
	std::vector< StatData > getMean( size_t id );

	//	Calculate the variance of data for a given event
	//	Variance - average of the squared differences between data points and the mean
	std::vector< StatData >	getVariance( size_t id );

	//	Sqrt of variance, also in units of the original data
	std::vector< StatData >	getStdDev( size_t id );

	/**
	 * \fn double getAverageTime(size_t id) const
	 * \return Return the arithmetic mean of all the samples that have been saved
	 */
	std::vector< StatData > getAverageTime( size_t id );

	/**
	 * \fn double getMinimumTime(size_t id) const
	 * \return Return the arithmetic min of all the samples that have been saved
	 */
	std::vector< StatData > getMinimumTime( size_t id );

	void queryOpenCL( size_t id );

	void ReleaseEvents();

public:
	/**
	 * \fn getInstance()
	 * \brief This returns a reference to the singleton timer.  Guarantees only 1 timer class is ever
	 *	instantiated within a compilable executable.
	 */
	static GpuStatTimer& getInstance( );

	/**
	 * \fn void Start( size_t id )
	 * \brief Start the timer
	 * \sa Stop(), Reset()
	 */
	void Start( size_t id );

	/**
	 * \fn void Stop( size_t id )
	 * \brief Stop the timer
	 * \sa Start(), Reset()
	 */
	void Stop( size_t id );

	/**
	 * \fn void AddSample( const cl_event ev )
	 * \brief Explicitely add a timing sample into the class
	 */
	virtual void AddSample( clfftPlanHandle plHandle, FFTPlan* plan, cl_kernel kern, cl_uint numQueuesAndEvents, cl_event* ev,
		const std::vector< size_t >& gWorkSize, const std::vector< size_t >& lWorkSize );

	/**
	 * \fn void Reset(void)
	 * \brief Reset the timer to 0
	 * \sa Start(), Stop()
	 */
	void Clear( );

	/**
	 * \fn void Reset(void)
	 * \brief Reset the timer to 0
	 * \sa Start(), Stop()
	 */
	void Reset( );

	void Reserve( size_t nEvents, size_t nSamples );

	size_t getUniqueID( const std::string& label, cl_uint groupID );

	//	Calculate the average/mean of data for a given event
	void	setNormalize( bool norm );

	void Print( );

	//	Using the stdDev of the entire population (of an id), eliminate those samples that fall
	//	outside some specified multiple of the stdDev.  This assumes that the population
	//	form a gaussian curve.
	size_t	pruneOutliers( cl_double multiple );
	std::vector< size_t > pruneOutliers( size_t id , cl_double multiple );
};

#endif // _STATISTICALTIMER_GPU_H_
