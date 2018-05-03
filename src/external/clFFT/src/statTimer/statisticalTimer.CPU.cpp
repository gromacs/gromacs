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


// StatTimer.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>
#include <functional>
#include <cmath>
#include <limits>

#include "statisticalTimer.CPU.h"
#include "../library/private.h"

#if defined( __GNUC__ )
	#include <sys/time.h>
#endif

//	Format an unsigned number with comma thousands separator
//
template< typename T >		// T could be 32-bit or 64-bit
std::basic_string<TCHAR> commatize (T number)
{
	static TCHAR scratch [8*sizeof(T)];

	register TCHAR * ptr = scratch + countOf( scratch );
	*(--ptr) = 0;

	for (int digits = 3; ; )
	{
		*(--ptr) = '0' + int (number % 10);
		number /= 10;
		if (0 == number)
			break;
		if (--digits <= 0)
		{
			*(--ptr) = ',';
			digits = 3;
		}
	}

	return std::basic_string<TCHAR> (ptr);
}

//	Functor object to help with accumulating values in vectors
template< typename T >
struct Accumulator: public std::unary_function< T, void >
{
	T acc;

	Accumulator( ): acc( 0 ) {}
	void operator( )(T x) { acc += x; }
};

//	Unary predicate used for remove_if() algorithm
//	Currently, RangeType is expected to be a floating point type, and ValType an integer type
template< typename RangeType, typename ValType >
struct PruneRange
{
	RangeType lower, upper;

	PruneRange( RangeType mean, RangeType stdev ): lower( mean-stdev ), upper( mean+stdev ) {}

	bool operator( )( ValType val )
	{
		//	These comparisons can be susceptible to signed/unsigned casting problems
		//	This is why we cast ValType to RangeType, because RangeType should always be floating and signed
		if( static_cast< RangeType >( val ) < lower )
			return true;
		else if( static_cast< RangeType >( val ) > upper )
			return true;

		return false;
	}
};

CpuStatTimer&
CpuStatTimer::getInstance( )
{
	static	CpuStatTimer	timer;
	return	timer;
}

CpuStatTimer::CpuStatTimer( ): nEvents( 0 ), nSamples( 0 ), normalize( true )
{
#if defined( _WIN32 )
	//	OS call to get ticks per second2
	::QueryPerformanceFrequency( reinterpret_cast<LARGE_INTEGER*>( &clkFrequency ) );
#else
	res.tv_sec	= 0;
	res.tv_nsec	= 0;
	clkFrequency 	= 0;

	//	clock_getres() return 0 for success
	//	If the function fails (monotonic clock not supported), we default to a lower resolution timer
//	if( ::clock_getres( CLOCK_MONOTONIC, &res ) )
	{
		clkFrequency = 1000000;
	}
//	else
//	{
//	    // Turn time into frequency
//		clkFrequency = res.tv_nsec * 1000000000;
//	}

#endif
}

CpuStatTimer::~CpuStatTimer( )
{}

void
CpuStatTimer::Clear( )
{
	labelID.clear( );
	clkStart.clear( );
	clkTicks.clear( );
}

void
CpuStatTimer::Reset( )
{
	if( nEvents == 0 || nSamples == 0 )
		throw	std::runtime_error( "StatisticalTimer::Reserve( ) was not called before Reset( )" );

	clkStart.clear( );
	clkTicks.clear( );

	clkStart.resize( nEvents );
	clkTicks.resize( nEvents );

	for( cl_uint	i = 0; i < nEvents; ++i )
	{
		clkTicks.at( i ).reserve( nSamples );
	}

	return;
}

//	The caller can pre-allocate memory, to improve performance.
//	nEvents is an approximate value for how many seperate events the caller will think
//	they will need, and nSamples is a hint on how many samples we think we will take
//	per event
void
CpuStatTimer::Reserve( size_t nEvents, size_t nSamples )
{
	this->nEvents	= std::max< size_t >( 1, nEvents );
	this->nSamples	= std::max< size_t >( 1, nSamples );

	Clear( );
	labelID.reserve( nEvents );

	clkStart.resize( nEvents );
	clkTicks.resize( nEvents );

	for( cl_uint i = 0; i < nEvents; ++i )
	{
		clkTicks.at( i ).reserve( nSamples );
	}
}

void
CpuStatTimer::setNormalize( bool norm )
{
	normalize = norm;
}

void
CpuStatTimer::Start( size_t id )
{
#if defined( _WIN32 )
	::QueryPerformanceCounter( reinterpret_cast<LARGE_INTEGER*>( &clkStart.at( id ) ) );
#else
	if( clkFrequency )
	{
		struct timeval s;
		gettimeofday( &s, 0 );
		clkStart.at( id ) = (cl_ulong)s.tv_sec * 1000000 + (cl_ulong)s.tv_usec;
	}
	else
	{

	}
#endif
}

void
CpuStatTimer::Stop( size_t id )
{
	cl_ulong n;

#if defined( _WIN32 )
	::QueryPerformanceCounter( reinterpret_cast<LARGE_INTEGER*>( &n ) );
#else
	struct timeval s;
	gettimeofday( &s, 0 );
	n = (cl_ulong)s.tv_sec * 1000000 + (cl_ulong)s.tv_usec;
#endif

	n		-= clkStart.at( id );
	clkStart.at( id )	= 0;
	AddSample( id, n );
}

void
CpuStatTimer::AddSample( const size_t id, const cl_ulong n )
{
	clkTicks.at( id ).push_back( n );
}

//	This function's purpose is to provide a mapping from a 'friendly' human readable text string
//	to an index into internal data structures.
size_t
CpuStatTimer::getUniqueID( const std::string& label, cl_uint groupID )
{
	//	I expect labelID will hardly ever grow beyond 30, so it's not of any use
	//	to keep this sorted and do a binary search

	labelPair	sItem	= std::make_pair( label, groupID );

	stringVector::iterator	iter;
	iter	= std::find( labelID.begin(), labelID.end(), sItem );

	if( iter != labelID.end( ) )
		return	std::distance( labelID.begin( ), iter );

	labelID.push_back( sItem );

	return	labelID.size( ) - 1;

}

cl_double
CpuStatTimer::getMean( size_t id ) const
{
	if( clkTicks.empty( ) )
		return	0;

	size_t	N	= clkTicks.at( id ).size( );

	Accumulator<cl_ulong> sum = std::for_each( clkTicks.at( id ).begin(), clkTicks.at( id ).end(), Accumulator<cl_ulong>() );

	return	static_cast<cl_double>( sum.acc ) / N;
}

cl_double
CpuStatTimer::getVariance( size_t id ) const
{
	if( clkTicks.empty( ) )
		return	0;

	cl_double	mean	= getMean( id );

	size_t	N	= clkTicks.at( id ).size( );
	cl_double	sum	= 0;

	for( cl_uint i = 0; i < N; ++i )
	{
		cl_double	diff	= clkTicks.at( id ).at( i ) - mean;
		diff	*= diff;
		sum		+= diff;
	}

	return	 sum / N;
}

cl_double
CpuStatTimer::getStdDev( size_t id ) const
{
	cl_double	variance	= getVariance( id );

	return	sqrt( variance );
}

cl_double
CpuStatTimer::getAverageTime( size_t id ) const
{
	if( normalize )
		return getMean( id ) / clkFrequency;
	else
		return getMean( id );
}

cl_double
CpuStatTimer::getMinimumTime( size_t id ) const
{
	clkVector::const_iterator iter	= std::min_element( clkTicks.at( id ).begin( ), clkTicks.at( id ).end( ) );

	if( iter != clkTicks.at( id ).end( ) )
	{
		if( normalize )
			return static_cast<cl_double>( *iter ) / clkFrequency;
		else
			return static_cast<cl_double>( *iter );
	}
	else
		return	0;
}

std::vector< size_t >
CpuStatTimer::pruneOutliers( size_t id , cl_double multiple )
{
	//if( clkTicks.empty( ) )
	//	return	std::vector< size_t >( );

	//cl_double	mean	= getMean( id );
	//cl_double	stdDev	= getStdDev( id );

	//clkVector&	clks = clkTicks.at( id );

	////	Look on p. 379, "The C++ Standard Library"
	////	std::remove_if does not actually erase, it only copies elements, it returns new 'logical' end
	//clkVector::iterator	newEnd	= std::remove_if( clks.begin( ), clks.end( ), PruneRange< cl_double,cl_ulong >( mean, multiple*stdDev ) );

	//clkVector::difference_type dist	= std::distance( newEnd, clks.end( ) );

	//if( dist != 0 )
	//	clks.erase( newEnd, clks.end( ) );

	//assert( dist < std::numeric_limits< cl_uint >::max( ) );

	return std::vector< size_t >( );
}

size_t
CpuStatTimer::pruneOutliers( cl_double multiple )
{
	size_t	tCount	= 0;

	//for( cl_uint l = 0; l < labelID.size( ); ++l )
	//{
	//	size_t lCount	= pruneOutliers( l , multiple );
	//	std::clog << "\tStatisticalTimer:: Pruning " << lCount << " samples from " << labelID[l].first << std::endl;
	//	tCount += lCount;
	//}

	return	tCount;
}

void
CpuStatTimer::Print( )
{
	//double flops = fFunc( );

	//for( cl_uint i = 0; i < labelID.size( ); ++i )
	//{
	//	double timeNs	= getAverageTime( i );
	//	double gFlops	= flops / timeNs;

	//	std::cout << labelID[ i ].first << std::endl;
	//	tout << std::setw( 10 ) << "Time:" << std::setw( 10 )  << commatize( static_cast< cl_ulong >( timeNs ) )
	//		<< " ns" << std::endl;
	//	tout << std::setw( 10 ) << "Gflops:" << std::setw( 10 ) << gFlops << std::endl;
	//}
}

//	Defining an output print operator
std::ostream&
operator<<( std::ostream& os, const CpuStatTimer& st )
{
	if( st.clkTicks.empty( ) )
		return	os;

	std::ios::fmtflags bckup	= os.flags( );

	for( cl_uint l = 0; l < st.labelID.size( ); ++l )
	{
		cl_ulong min	= 0;
		CpuStatTimer::clkVector::const_iterator iter	= std::min_element( st.clkTicks.at( l ).begin( ), st.clkTicks.at( l ).end( ) );

		if( iter != st.clkTicks.at( l ).end( ) )
			min		= *iter;

		os << st.labelID[l].first << ", " << st.labelID[l].second << std::fixed << std::endl;
		os << "Min:," << min << std::endl;
		os << "Mean:," << st.getMean( l ) << std::endl;
		os << "StdDev:," << st.getStdDev( l ) << std::endl;
		os << "AvgTime:," << st.getAverageTime( l ) << std::endl;
		os << "MinTime:," << st.getMinimumTime( l ) << std::endl;

		//for( cl_uint	t = 0; t < st.clkTicks[l].size( ); ++t )
		//{
		//	os << st.clkTicks[l][t]<< ",";
		//}
		os << "\n" << std::endl;

	}

	os.flags( bckup );

	return	os;
}
