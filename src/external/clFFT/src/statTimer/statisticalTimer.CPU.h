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
#ifndef _STATISTICALTIMER_CPU_H_
#define _STATISTICALTIMER_CPU_H_
#include <iosfwd>
#include <vector>
#include <algorithm>
#ifdef __FreeBSD__
#include <sys/timespec.h>
#endif
#include "statisticalTimer.h"

/**
 * \file clfft.StatisticalTimer.CPU.h
 * \brief A timer class that provides a cross platform timer for use
 * in timing code progress with a high degree of accuracy.
 *	This class is implemented entirely in the header, to facilitate inclusion into multiple
 *	projects without needing to compile an object file for each project.
 */

class CpuStatTimer : public baseStatTimer
{
	//	Private typedefs
	typedef std::vector< cl_ulong > clkVector;
	typedef	std::pair< std::string, cl_uint > labelPair;
	typedef	std::vector< labelPair > stringVector;

	//	In order to calculate statistics <std. dev.>, we need to keep a history of our timings
	stringVector	labelID;
	clkVector	clkStart;
	std::vector< clkVector >	clkTicks;

	//	How many clockticks in a second
	cl_ulong	clkFrequency;

	//	For linux; the resolution of a high-precision timer
	//  Mingw32 does not define timespec; can use windows timers
#if !defined( _WIN32 )
	timespec res;
#endif

	//	Saved sizes for our vectors, used in Reset() to reallocate vectors
	clkVector::size_type	nEvents, nSamples;

	//	This setting controls whether the Timer should convert samples into time by dividing by the
	//	clock frequency
	bool normalize;

	/**
	 * \fn StatisticalTimer()
	 * \brief Constructor for StatisticalTimer that initializes the class
	 *	This is private so that user code cannot create their own instantiation.  Instead, you
	 *	must go through getInstance( ) to get a reference to the class.
	 */
	CpuStatTimer( );

	/**
	 * \fn ~StatisticalTimer()
	 * \brief Destructor for StatisticalTimer that cleans up the class
	 */
	~CpuStatTimer( );

	/**
	 * \fn StatisticalTimer(const StatisticalTimer& )
	 * \brief Copy constructors do not make sense for a singleton, disallow copies
	 */
	CpuStatTimer( const CpuStatTimer& );

	/**
	 * \fn operator=( const StatisticalTimer& )
	 * \brief Assignment operator does not make sense for a singleton, disallow assignments
	 */
	CpuStatTimer& operator=( const CpuStatTimer& );

	friend std::ostream& operator<<( std::ostream& os, const CpuStatTimer& s );

	/**
	 * \fn void AddSample( const size_t id, const cl_ulong n )
	 * \brief Explicitely add a timing sample into the class
	 */
	void AddSample( const size_t id, const cl_ulong n );

	//	Calculate the average/mean of data for a given event
	cl_double	getMean( size_t id ) const;

	//	Calculate the variance of data for a given event
	//	Variance - average of the squared differences between data points and the mean
	cl_double	getVariance( size_t id ) const;

	//	Sqrt of variance, also in units of the original data
	cl_double	getStdDev( size_t id ) const;

	/**
	 * \fn double getAverageTime(size_t id) const
	 * \return Return the arithmetic mean of all the samples that have been saved
	 */
	cl_double getAverageTime( size_t id ) const;

	/**
	 * \fn double getMinimumTime(size_t id) const
	 * \return Return the arithmetic min of all the samples that have been saved
	 */
	cl_double getMinimumTime( size_t id ) const;

public:
	/**
	 * \fn getInstance()
	 * \brief This returns a reference to the singleton timer.  Guarantees only 1 timer class is ever
	 *	instantiated within a compilable executable.
	 */
	static CpuStatTimer& getInstance( );

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

#endif // _STATISTICALTIMER_CPU_H_
