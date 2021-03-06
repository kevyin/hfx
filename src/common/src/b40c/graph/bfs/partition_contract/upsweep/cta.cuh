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
 * Thanks!
 * 
 ******************************************************************************/


/******************************************************************************
 * CTA tile-processing abstraction for multi-GPU BFS frontier
 * contraction+binning (upsweep reduction)
 ******************************************************************************/

#pragma once

#include <b40c/partition/upsweep/cta.cuh>

#include <b40c/graph/bfs/partition_contract/upsweep/tile.cuh>

namespace b40c {
namespace graph {
namespace bfs {
namespace partition_contract {
namespace upsweep {


/**
 * CTA tile-processing abstraction for multi-GPU BFS frontier
 * contraction+binning (upsweep reduction)
 *
 * Derives from partition::upsweep::Cta
 */
template <typename KernelPolicy>
struct Cta :
	partition::upsweep::Cta<
		KernelPolicy,
		Cta<KernelPolicy>,			// This class
		Tile>						// radix_sort::upsweep::Tile
{
	//---------------------------------------------------------------------
	// Typedefs and Constants
	//---------------------------------------------------------------------

	// Base class type
	typedef partition::upsweep::Cta<KernelPolicy, Cta, Tile> Base;

	typedef typename KernelPolicy::SmemStorage 				SmemStorage;
	typedef typename KernelPolicy::VertexId 				VertexId;
	typedef typename KernelPolicy::ValidFlag				ValidFlag;
	typedef typename KernelPolicy::VisitedMask 			VisitedMask;

	typedef typename KernelPolicy::KeyType 					KeyType;
	typedef typename KernelPolicy::SizeT 					SizeT;


	//---------------------------------------------------------------------
	// Members
	//---------------------------------------------------------------------

	// Input and output device pointers
	ValidFlag				*d_flags_out;
	VisitedMask 			*d_visited_mask;

	// Smem storage for reduction tree and hashing scratch
	volatile VertexId 		(*vid_hashtable)[SmemStorage::WARP_HASH_ELEMENTS];
	volatile VertexId		*history;

	// Number of GPUs to partition the frontier into
	int num_gpus;


	//---------------------------------------------------------------------
	// Methods
	//---------------------------------------------------------------------

	/**
	 * Constructor
	 */
	__device__ __forceinline__ Cta(
		SmemStorage 	&smem_storage,
		int 			num_gpus,
		VertexId 		*d_in,
		ValidFlag		*d_flags_out,
		SizeT 			*d_spine,
		VisitedMask 	*d_visited_mask) :
			Base(
				smem_storage,
				d_in, 				// d_in_keys
				d_spine),
			d_flags_out(d_flags_out),
			d_visited_mask(d_visited_mask),
			vid_hashtable(smem_storage.vid_hashtable),
			history(smem_storage.history),
			num_gpus(num_gpus)
	{
		// Initialize history filter
		for (int offset = threadIdx.x; offset < SmemStorage::HISTORY_HASH_ELEMENTS; offset += KernelPolicy::THREADS) {
			history[offset] = -1;
		}
	}

};



} // namespace upsweep
} // namespace partition_contract
} // namespace bfs
} // namespace graph
} // namespace b40c

