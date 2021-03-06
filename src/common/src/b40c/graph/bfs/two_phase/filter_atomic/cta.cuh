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
 * CTA tile-processing abstraction for BFS edge-frontier filtering
 ******************************************************************************/

#pragma once

#include <b40c/util/io/modified_load.cuh>
#include <b40c/util/io/modified_store.cuh>
#include <b40c/util/io/initialize_tile.cuh>
#include <b40c/util/io/load_tile.cuh>
#include <b40c/util/io/store_tile.cuh>
#include <b40c/util/io/scatter_tile.cuh>
#include <b40c/util/cta_work_distribution.cuh>

#include <b40c/util/operators.cuh>
#include <b40c/util/reduction/serial_reduce.cuh>
#include <b40c/util/reduction/tree_reduce.cuh>

namespace b40c {
namespace graph {
namespace bfs {
namespace two_phase {
namespace filter_atomic {


/**
 * CTA tile-processing abstraction for BFS edge-frontier filtering
 */
template <typename KernelPolicy>
struct Cta
{
	//---------------------------------------------------------------------
	// Typedefs and Constants
	//---------------------------------------------------------------------

	typedef typename KernelPolicy::VertexId 		VertexId;
	typedef typename KernelPolicy::ValidFlag		ValidFlag;
	typedef typename KernelPolicy::VisitedMask 		VisitedMask;
	typedef typename KernelPolicy::SizeT 			SizeT;

	typedef typename KernelPolicy::RakingDetails 	RakingDetails;
	typedef typename KernelPolicy::SmemStorage		SmemStorage;

	//---------------------------------------------------------------------
	// Members
	//---------------------------------------------------------------------

	// Input and output device pointers
	VertexId 				*d_in;						// Incoming edge frontier
	VertexId 				*d_out;						// Outgoing filtered edge frontier
	VertexId 				*d_predecessor_in;			// Incoming predecessor edge frontier (used when KernelPolicy::MARK_PREDECESSORS)
	VertexId 				*d_predecessor_out;			// Outgoing predecessor edge frontier (used when KernelPolicy::MARK_PREDECESSORS)
	VisitedMask 			*d_visited_mask;			// Mask for detecting visited status

	// Work progress
	VertexId 				queue_index;				// Current frontier queue counter index
	util::CtaWorkProgress	&work_progress;				// Atomic workstealing and queueing counters
	SizeT					max_filtered_frontier;		// Maximum size (in elements) of outgoing filtered edge frontier

	// Operational details for raking scan grid
	RakingDetails 			raking_details;

	// Shared memory for the CTA
	SmemStorage				&smem_storage;



	//---------------------------------------------------------------------
	// Helper Structures
	//---------------------------------------------------------------------

	/**
	 * Tile of incoming edge frontier to process
	 */
	template <
		int LOG_LOADS_PER_TILE,
		int LOG_LOAD_VEC_SIZE>
	struct Tile
	{
		//---------------------------------------------------------------------
		// Typedefs and Constants
		//---------------------------------------------------------------------

		enum {
			LOADS_PER_TILE 		= 1 << LOG_LOADS_PER_TILE,
			LOAD_VEC_SIZE 		= 1 << LOG_LOAD_VEC_SIZE
		};


		//---------------------------------------------------------------------
		// Members
		//---------------------------------------------------------------------

		// Dequeued vertex ids
		VertexId 	vertex_id[LOADS_PER_TILE][LOAD_VEC_SIZE];
		VertexId 	predecessor_id[LOADS_PER_TILE][LOAD_VEC_SIZE];

		// Whether or not the corresponding vertex_id is valid for exploring
		ValidFlag 	flags[LOADS_PER_TILE][LOAD_VEC_SIZE];

		// Global scatter offsets
		SizeT 		ranks[LOADS_PER_TILE][LOAD_VEC_SIZE];

		//---------------------------------------------------------------------
		// Helper Structures
		//---------------------------------------------------------------------

		/**
		 * Iterate over vertex id
		 */
		template <int LOAD, int VEC, int dummy = 0>
		struct Iterate
		{
			/**
			 * InitFlags
			 */
			static __device__ __forceinline__ void InitFlags(Tile *tile)
			{
				// Initially valid if vertex-id is valid
				tile->flags[LOAD][VEC] = (tile->vertex_id[LOAD][VEC] == -1) ? 0 : 1;
				tile->ranks[LOAD][VEC] = tile->flags[LOAD][VEC];

				// Next
				Iterate<LOAD, VEC + 1>::InitFlags(tile);
			}

			/**
			 * BitmaskCull
			 */
			static __device__ __forceinline__ void BitmaskCull(
				Cta *cta,
				Tile *tile)
			{
				if (tile->vertex_id[LOAD][VEC] != -1) {

					// Location of mask byte to read
					SizeT mask_byte_offset = (tile->vertex_id[LOAD][VEC] & KernelPolicy::VERTEX_ID_MASK) >> 3;

					// Bit in mask byte corresponding to current vertex id
					VisitedMask mask_bit = 1 << (tile->vertex_id[LOAD][VEC] & 7);

					// Read byte from from visited mask in tex
					VisitedMask tex_mask_byte = tex1Dfetch(
						contract_atomic::BitmaskTex<VisitedMask>::ref,
						mask_byte_offset);

					if (mask_bit & tex_mask_byte) {

						// Seen it
						tile->vertex_id[LOAD][VEC] = -1;

					} else {

						VisitedMask mask_byte;
						util::io::ModifiedLoad<util::io::ld::cg>::Ld(
							mask_byte, cta->d_visited_mask + mask_byte_offset);

						mask_byte |= tex_mask_byte;

						if (mask_bit & mask_byte) {

							// Seen it
							tile->vertex_id[LOAD][VEC] = -1;

						} else {

							// Update with best effort
							mask_byte |= mask_bit;
							util::io::ModifiedStore<util::io::st::cg>::St(
								mask_byte,
								cta->d_visited_mask + mask_byte_offset);
						}
					}
				}

				// Next
				Iterate<LOAD, VEC + 1>::BitmaskCull(cta, tile);
			}
		};


		/**
		 * Iterate next load
		 */
		template <int LOAD, int dummy>
		struct Iterate<LOAD, LOAD_VEC_SIZE, dummy>
		{
			// InitFlags
			static __device__ __forceinline__ void InitFlags(Tile *tile)
			{
				Iterate<LOAD + 1, 0>::InitFlags(tile);
			}

			// BitmaskCull
			static __device__ __forceinline__ void BitmaskCull(Cta *cta, Tile *tile)
			{
				Iterate<LOAD + 1, 0>::BitmaskCull(cta, tile);
			}
		};


		/**
		 * Terminate iteration
		 */
		template <int dummy>
		struct Iterate<LOADS_PER_TILE, 0, dummy>
		{
			// InitFlags
			static __device__ __forceinline__ void InitFlags(Tile *tile) {}

			// BitmaskCull
			static __device__ __forceinline__ void BitmaskCull(Cta *cta, Tile *tile) {}
		};


		//---------------------------------------------------------------------
		// Interface
		//---------------------------------------------------------------------

		/**
		 * Initializer
		 */
		__device__ __forceinline__ void InitFlags()
		{
			Iterate<0, 0>::InitFlags(this);
		}

		/**
		 * Culls vertices based upon whether or not we've set a bit for them
		 * in the d_visited_mask bitmask
		 */
		__device__ __forceinline__ void BitmaskCull(Cta *cta)
		{
			Iterate<0, 0>::BitmaskCull(cta, this);
		}

	};




	//---------------------------------------------------------------------
	// Methods
	//---------------------------------------------------------------------

	/**
	 * Constructor
	 */
	__device__ __forceinline__ Cta(
		VertexId 				queue_index,
		SmemStorage 			&smem_storage,
		VertexId 				*d_in,
		VertexId 				*d_out,
		VertexId 				*d_predecessor_in,
		VertexId 				*d_predecessor_out,
		VisitedMask 			*d_visited_mask,
		util::CtaWorkProgress	&work_progress,
		SizeT					max_filtered_frontier) :

			queue_index(queue_index),
			raking_details(
				smem_storage.state.raking_elements,
				smem_storage.state.warpscan,
				0),
			smem_storage(smem_storage),
			d_in(d_in),
			d_out(d_out),
			d_predecessor_in(d_predecessor_in),
			d_predecessor_out(d_predecessor_out),
			d_visited_mask(d_visited_mask),
			work_progress(work_progress),
			max_filtered_frontier(max_filtered_frontier)
	{
	}


	/**
	 * Process a single, full tile
	 */
	__device__ __forceinline__ void ProcessTile(
		SizeT cta_offset,
		const SizeT &guarded_elements = KernelPolicy::TILE_ELEMENTS)
	{
		Tile<KernelPolicy::LOG_LOADS_PER_TILE, KernelPolicy::LOG_LOAD_VEC_SIZE> tile;

		// Load tile
		util::io::LoadTile<
			KernelPolicy::LOG_LOADS_PER_TILE,
			KernelPolicy::LOG_LOAD_VEC_SIZE,
			KernelPolicy::THREADS,
			KernelPolicy::QUEUE_READ_MODIFIER,
			false>::LoadValid(
				tile.vertex_id,
				d_in,
				cta_offset,
				guarded_elements,
				(VertexId) -1);

		if (KernelPolicy::MARK_PREDECESSORS) {

			// Load predecessor vertices as well
			util::io::LoadTile<
				KernelPolicy::LOG_LOADS_PER_TILE,
				KernelPolicy::LOG_LOAD_VEC_SIZE,
				KernelPolicy::THREADS,
				KernelPolicy::QUEUE_READ_MODIFIER,
				false>::LoadValid(
					tile.predecessor_id,
					d_predecessor_in,
					cta_offset,
					guarded_elements);
		}

		// Cull visited vertices and update discovered vertices using global visited mask
		tile.BitmaskCull(this);

		// Init valid flags and ranks
		tile.InitFlags();

		// Scan tile of ranks, using an atomic add to reserve
		// space in the contracted queue, seeding ranks
		util::Sum<SizeT> scan_op;
		SizeT new_queue_offset = util::scan::CooperativeTileScan<KernelPolicy::LOAD_VEC_SIZE>::ScanTileWithEnqueue(
			raking_details,
			tile.ranks,
			work_progress.GetQueueCounter<SizeT>(queue_index + 1),
			scan_op);

		// Check updated queue offset for overflow due to redundant expansion
		if (new_queue_offset >= max_filtered_frontier) {
			work_progress.SetOverflow<SizeT>();
			util::ThreadExit();
		}

		// Scatter directly (without first contracting in smem scratch), predicated
		// on flags
		util::io::ScatterTile<
			KernelPolicy::LOG_LOADS_PER_TILE,
			KernelPolicy::LOG_LOAD_VEC_SIZE,
			KernelPolicy::THREADS,
			KernelPolicy::QUEUE_WRITE_MODIFIER>::Scatter(
				d_out,
				tile.vertex_id,
				tile.flags,
				tile.ranks);

		if (KernelPolicy::MARK_PREDECESSORS) {

			util::io::ScatterTile<
				KernelPolicy::LOG_LOADS_PER_TILE,
				KernelPolicy::LOG_LOAD_VEC_SIZE,
				KernelPolicy::THREADS,
				KernelPolicy::QUEUE_WRITE_MODIFIER>::Scatter(
					d_predecessor_out,
					tile.predecessor_id,
					tile.flags,
					tile.ranks);
		}

	}
};


} // namespace filter_atomic
} // namespace two_phase
} // namespace bfs
} // namespace graph
} // namespace b40c

