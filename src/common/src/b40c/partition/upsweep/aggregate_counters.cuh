/******************************************************************************
 * 
 * Copyright 2010-2012 Duane Merrill
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. 
 *
 * For more information, see our Google Code project site: 
 * http://code.google.com/p/back40computing/
 * 
 ******************************************************************************/

/******************************************************************************
 * Aggregate-counter functionality for partitioning upsweep reduction kernels
 ******************************************************************************/

#pragma once

namespace b40c {
namespace partition {
namespace upsweep {


/**
 * Thread-local aggregate counters.  Each warp must periodically expand and
 * aggregate its shared-memory composite-counter lanes into registers before
 * they overflow.
 *
 * Each thread will aggregate a segment of (lane-length / warp-size)
 * composite-counters within lane belonging to its warp.  There
 * are four encoded bins within each composite counter.
 */
template <typename KernelPolicy>
struct AggregateCounters
{
	//---------------------------------------------------------------------
	// Typedefs and constants
	//---------------------------------------------------------------------

	typedef typename KernelPolicy::SizeT SizeT;

	enum {
		LANES_PER_WARP 						= KernelPolicy::LANES_PER_WARP,
		COMPOSITES_PER_LANE_PER_THREAD 		= KernelPolicy::COMPOSITES_PER_LANE_PER_THREAD,
		WARPS								= KernelPolicy::WARPS,
		COMPOSITE_LANES						= KernelPolicy::COMPOSITE_LANES,
	};

	//---------------------------------------------------------------------
	// Members
	//---------------------------------------------------------------------

	// Aggregate counters (four encoded bins per lane by the number of lanes
	// aggregated per warp)
	SizeT local_counts[KernelPolicy::LANES_PER_WARP][4];


	//---------------------------------------------------------------------
	// Iteration Structures
	//---------------------------------------------------------------------

	/**
	 * Iterate next composite counter
	 */
	template <int WARP_LANE, int THREAD_COMPOSITE, int dummy = 0>
	struct Iterate
	{
		// ExtractComposites
		template <typename Cta, typename AggregateCounters>
		static __device__ __forceinline__ void ExtractComposites(
			Cta *cta,
			AggregateCounters *aggregate_counters)
		{
			int lane				= (WARP_LANE * WARPS) + cta->warp_id;
			int composite			= (THREAD_COMPOSITE * B40C_WARP_THREADS(__B40C_CUDA_ARCH__)) + cta->warp_idx;

			aggregate_counters->local_counts[WARP_LANE][0] += cta->smem_storage.composite_counters.counters[lane][composite][0];
			aggregate_counters->local_counts[WARP_LANE][1] += cta->smem_storage.composite_counters.counters[lane][composite][1];
			aggregate_counters->local_counts[WARP_LANE][2] += cta->smem_storage.composite_counters.counters[lane][composite][2];
			aggregate_counters->local_counts[WARP_LANE][3] += cta->smem_storage.composite_counters.counters[lane][composite][3];

			Iterate<WARP_LANE, THREAD_COMPOSITE + 1>::ExtractComposites(cta, aggregate_counters);
		}
	};

	/**
	 * Iterate next lane
	 */
	template <int WARP_LANE, int dummy>
	struct Iterate<WARP_LANE, COMPOSITES_PER_LANE_PER_THREAD, dummy>
	{
		// ExtractComposites
		template <typename Cta, typename AggregateCounters>
		static __device__ __forceinline__ void ExtractComposites(
			Cta *cta,
			AggregateCounters *aggregate_counters)
		{
			Iterate<WARP_LANE + 1, 0>::ExtractComposites(cta, aggregate_counters);
		}

		// ShareCounters
		template <typename Cta, typename AggregateCounters>
		static __device__ __forceinline__ void ShareCounters(
			Cta *cta,
			AggregateCounters *aggregate_counters)
		{
			int lane				= (WARP_LANE * WARPS) + cta->warp_id;
			int row 				= lane << 2;	// lane * 4;

			cta->smem_storage.aggregate[row + 0][cta->warp_idx] = aggregate_counters->local_counts[WARP_LANE][0];
			cta->smem_storage.aggregate[row + 1][cta->warp_idx] = aggregate_counters->local_counts[WARP_LANE][1];
			cta->smem_storage.aggregate[row + 2][cta->warp_idx] = aggregate_counters->local_counts[WARP_LANE][2];
			cta->smem_storage.aggregate[row + 3][cta->warp_idx] = aggregate_counters->local_counts[WARP_LANE][3];

			Iterate<WARP_LANE + 1, COMPOSITES_PER_LANE_PER_THREAD>::ShareCounters(cta, aggregate_counters);
		}

		// ResetCounters
		template <typename Cta, typename AggregateCounters>
		static __device__ __forceinline__ void ResetCounters(
			Cta *cta,
			AggregateCounters *aggregate_counters)
		{
			aggregate_counters->local_counts[WARP_LANE][0] = 0;
			aggregate_counters->local_counts[WARP_LANE][1] = 0;
			aggregate_counters->local_counts[WARP_LANE][2] = 0;
			aggregate_counters->local_counts[WARP_LANE][3] = 0;

			Iterate<WARP_LANE + 1, COMPOSITES_PER_LANE_PER_THREAD>::ResetCounters(cta, aggregate_counters);
		}
	};

	/**
	 * Terminate iteration
	 */
	template <int dummy>
	struct Iterate<LANES_PER_WARP, 0, dummy>
	{
		// ExtractComposites
		template <typename Cta, typename AggregateCounters>
		static __device__ __forceinline__ void ExtractComposites(
			Cta *cta, AggregateCounters *aggregate_counters) {}
	};

	/**
	 * Terminate iteration
	 */
	template <int dummy>
	struct Iterate<LANES_PER_WARP, COMPOSITES_PER_LANE_PER_THREAD, dummy>
	{
		// ShareCounters
		template <typename Cta, typename AggregateCounters>
		static __device__ __forceinline__ void ShareCounters(
			Cta *cta, AggregateCounters *aggregate_counters) {}

		// ResetCounters
		template <typename Cta, typename AggregateCounters>
		static __device__ __forceinline__ void ResetCounters(
			Cta *cta, AggregateCounters *aggregate_counters) {}
	};

	//---------------------------------------------------------------------
	// Interface
	//---------------------------------------------------------------------

	/**
	 * Extracts and aggregates the shared-memory composite counters for each
	 * composite-counter lane owned by this warp
	 */
	template <typename Cta>
	__device__ __forceinline__ void ExtractComposites(Cta *cta)
	{
		if (cta->warp_id < COMPOSITE_LANES) {
			Iterate<0, 0>::ExtractComposites(cta, this);
		}
	}

	/**
	 * Places aggregate-counters into shared storage for final bin-wise reduction
	 */
	template <typename Cta>
	__device__ __forceinline__ void ShareCounters(Cta *cta)
	{
		if (cta->warp_id < COMPOSITE_LANES) {
			Iterate<0, COMPOSITES_PER_LANE_PER_THREAD>::ShareCounters(cta, this);
		}
	}

	/**
	 * Resets the aggregate counters
	 */
	template <typename Cta>
	__device__ __forceinline__ void ResetCounters(Cta *cta)
	{
		Iterate<0, COMPOSITES_PER_LANE_PER_THREAD>::ResetCounters(cta, this);
	}
};



} // namespace upsweep
} // namespace partition
} // namespace b40c

