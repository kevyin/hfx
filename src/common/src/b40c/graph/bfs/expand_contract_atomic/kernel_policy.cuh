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
 * Kernel configuration policy for BFS frontier expansion+contraction kernels.
 ******************************************************************************/

#pragma once

#include <b40c/util/basic_utils.cuh>
#include <b40c/util/cuda_properties.cuh>
#include <b40c/util/cta_work_distribution.cuh>
#include <b40c/util/soa_tuple.cuh>
#include <b40c/util/srts_grid.cuh>
#include <b40c/util/srts_details.cuh>
#include <b40c/util/io/modified_load.cuh>
#include <b40c/util/io/modified_store.cuh>

namespace b40c {
namespace graph {
namespace bfs {
namespace expand_contract_atomic {

/**
 * Kernel configuration policy for BFS frontier expansion+contraction kernels.
 *
 * Parameterizations of this type encapsulate our kernel-tuning parameters
 * (i.e., they are reflected via the static fields).
 *
 * Kernels can be specialized for problem-type, SM-version, etc. by parameterizing
 * them with different performance-tuned parameterizations of this type.  By
 * incorporating this type into the kernel code itself, we guide the compiler in
 * expanding/unrolling the kernel code for specific architectures and problem
 * types.
 */
template <
	// ProblemType type parameters
	typename _ProblemType,								// BFS problem type (e.g., b40c::graph::bfs::ProblemType)

	// Machine parameters
	int _CUDA_ARCH,										// CUDA SM architecture to generate code for

	// Behavioral control parameters
	bool _INSTRUMENT,									// Whether or not to record per-CTA clock timing statistics (for detecting load imbalance)
	int _SATURATION_QUIT,								// Early-exit threshold factor for the incoming frontier, i.e., quit when: frontier_size < SATURATION_QUIT * grid_size * TILE_SIZE.  (Can be -1 to never perform bitmask filtering)

	// Tunable parameters
	int _MIN_CTA_OCCUPANCY,												// Lower bound on number of CTAs to have resident per SM (influences per-CTA smem cache sizes and register allocation/spills)
	int _LOG_THREADS,													// Number of threads per CTA (log)
	int _LOG_LOAD_VEC_SIZE,												// Number of incoming frontier vertex-ids to dequeue in a single load (log)
	int _LOG_LOADS_PER_TILE,											// Number of such loads that constitute a tile of incoming frontier vertex-ids (log)
	int _LOG_RAKING_THREADS,											// Number of raking threads to use for prefix sum (log), range [5, LOG_THREADS]
	util::io::ld::CacheModifier _QUEUE_READ_MODIFIER,					// Load instruction cache-modifier for reading incoming frontier vertex-ids. Valid on SM2.0 or newer, where util::io::ld::cg is req'd for fused-iteration implementations incorporating software global barriers.
	util::io::ld::CacheModifier _COLUMN_READ_MODIFIER,					// Load instruction cache-modifier for reading CSR column-indices
	util::io::ld::CacheModifier _ROW_OFFSET_ALIGNED_READ_MODIFIER,		// Load instruction cache-modifier for reading CSR row-offsets (when 8-byte aligned)
	util::io::ld::CacheModifier _ROW_OFFSET_UNALIGNED_READ_MODIFIER,	// Load instruction cache-modifier for reading CSR row-offsets (when 4-byte aligned)
	util::io::st::CacheModifier _QUEUE_WRITE_MODIFIER,					// Store instruction cache-modifier for writing outgoign frontier vertex-ids. Valid on SM2.0 or newer, where util::io::st::cg is req'd for fused-iteration implementations incorporating software global barriers.
	bool _WORK_STEALING,												// Whether or not incoming frontier tiles are distributed via work-stealing or by even-share.
	int _WARP_GATHER_THRESHOLD,											// Adjacency-list length above which we expand an that list using coarser-grained warp-based cooperative expansion (below which we perform fine-grained scan-based expansion)
	int _CTA_GATHER_THRESHOLD,											// Adjacency-list length above which we expand an that list using coarsest-grained CTA-based cooperative expansion (below which we perform warp-based expansion)
	int _END_BITMASK_CULL,												// BFS Iteration after which to skip bitmask filtering (Alternatively 0 to never perform bitmask filtering, -1 to always perform bitmask filtering)
	int _LOG_SCHEDULE_GRANULARITY>										// The scheduling granularity of incoming frontier tiles (for even-share work distribution only) (log)

struct KernelPolicy : _ProblemType
{
	//---------------------------------------------------------------------
	// Constants and typedefs
	//---------------------------------------------------------------------

	typedef _ProblemType 					ProblemType;
	typedef typename ProblemType::VertexId 	VertexId;
	typedef typename ProblemType::SizeT 	SizeT;

	static const util::io::ld::CacheModifier QUEUE_READ_MODIFIER 					= _QUEUE_READ_MODIFIER;
	static const util::io::ld::CacheModifier COLUMN_READ_MODIFIER 					= _COLUMN_READ_MODIFIER;
	static const util::io::ld::CacheModifier ROW_OFFSET_ALIGNED_READ_MODIFIER 		= _ROW_OFFSET_ALIGNED_READ_MODIFIER;
	static const util::io::ld::CacheModifier ROW_OFFSET_UNALIGNED_READ_MODIFIER 	= _ROW_OFFSET_UNALIGNED_READ_MODIFIER;
	static const util::io::st::CacheModifier QUEUE_WRITE_MODIFIER 					= _QUEUE_WRITE_MODIFIER;

	enum {

		CUDA_ARCH 						= _CUDA_ARCH,
		INSTRUMENT						= _INSTRUMENT,
		SATURATION_QUIT					= _SATURATION_QUIT,

		LOG_THREADS 					= _LOG_THREADS,
		THREADS							= 1 << LOG_THREADS,

		LOG_LOAD_VEC_SIZE  				= _LOG_LOAD_VEC_SIZE,
		LOAD_VEC_SIZE					= 1 << LOG_LOAD_VEC_SIZE,

		LOG_LOADS_PER_TILE 				= _LOG_LOADS_PER_TILE,
		LOADS_PER_TILE					= 1 << LOG_LOADS_PER_TILE,

		LOG_LOAD_STRIDE					= LOG_THREADS + LOG_LOAD_VEC_SIZE,
		LOAD_STRIDE						= 1 << LOG_LOAD_STRIDE,

		LOG_RAKING_THREADS				= _LOG_RAKING_THREADS,
		RAKING_THREADS					= 1 << LOG_RAKING_THREADS,

		LOG_WARPS						= LOG_THREADS - B40C_LOG_WARP_THREADS(CUDA_ARCH),
		WARPS							= 1 << LOG_WARPS,

		LOG_TILE_ELEMENTS_PER_THREAD	= LOG_LOAD_VEC_SIZE + LOG_LOADS_PER_TILE,
		TILE_ELEMENTS_PER_THREAD		= 1 << LOG_TILE_ELEMENTS_PER_THREAD,

		LOG_TILE_ELEMENTS 				= LOG_TILE_ELEMENTS_PER_THREAD + LOG_THREADS,
		TILE_ELEMENTS					= 1 << LOG_TILE_ELEMENTS,

		LOG_SCHEDULE_GRANULARITY		= _LOG_SCHEDULE_GRANULARITY,
		SCHEDULE_GRANULARITY			= 1 << LOG_SCHEDULE_GRANULARITY,

		WORK_STEALING					= _WORK_STEALING,
		WARP_GATHER_THRESHOLD			= _WARP_GATHER_THRESHOLD,
		CTA_GATHER_THRESHOLD			= _CTA_GATHER_THRESHOLD,
		END_BITMASK_CULL 			= _END_BITMASK_CULL,
	};


	// Prefix sum raking grid for expansion allocations
	typedef util::RakingGrid<
		CUDA_ARCH,
		SizeT,									// Partial type (valid counts)
		LOG_THREADS,							// Depositing threads (the CTA size)
		LOG_LOADS_PER_TILE,						// Lanes (the number of loads)
		LOG_RAKING_THREADS,						// Raking threads
		true>									// There are prefix dependences between lanes
			RakingExpandGrid;


	// Operational details type for RakingExpandGrid raking grid type
	typedef util::RakingDetails<RakingExpandGrid> RakingExpandDetails;


	// Prefix sum raking grid for contraction allocations
	typedef util::RakingGrid<
		CUDA_ARCH,
		SizeT,									// Partial type (valid counts)
		LOG_THREADS,							// Depositing threads (the CTA size)
		0,										// 1 lane
		LOG_RAKING_THREADS,						// Raking threads
		true>									// There are prefix dependences between lanes
			RakingContractGrid;


	// Operational details type for RakingContractGrid raking grid type
	typedef util::RakingDetails<RakingContractGrid> RakingContractDetails;


	/**
	 * Shared memory storage type for the CTA
	 */
	struct SmemStorage
	{
		// Persistent shared state for the CTA
		struct State {
			// Type describing four shared memory channels per warp for intra-warp communication
			typedef SizeT 						WarpComm[WARPS][4];

			// Shared work-processing limits
			util::CtaWorkDistribution<SizeT>	work_decomposition;

			// Shared memory channels for intra-warp communication
			volatile WarpComm					warp_comm;
			int									cta_comm;

			// Storage for scanning local contract-expand ranks
			volatile SizeT 						warpscan[2][B40C_WARP_THREADS(CUDA_ARCH)];
			SizeT 								contract_raking_elements[RakingContractGrid::TOTAL_RAKING_ELEMENTS];

		} state;

		enum {
			// Amount of storage we can use for hashing scratch space under target occupancy
			MAX_SCRATCH_BYTES_PER_CTA		= (B40C_SMEM_BYTES(CUDA_ARCH) / _MIN_CTA_OCCUPANCY)
												- sizeof(State)
												- 128,											// Fudge-factor to guarantee occupancy

			SCRATCH_ELEMENT_SIZE 			= (ProblemType::MARK_PREDECESSORS) ?
													sizeof(SizeT) + sizeof(VertexId) :			// Need both gather offset and predecessor
													sizeof(SizeT),								// Just gather offset

			HASH_ELEMENTS					= MAX_SCRATCH_BYTES_PER_CTA / sizeof(VertexId),
			WARP_HASH_ELEMENTS				= 128,												// Collision hash table size (per warp)

			OFFSET_ELEMENTS					= MAX_SCRATCH_BYTES_PER_CTA / SCRATCH_ELEMENT_SIZE,
			PARENT_ELEMENTS					= (ProblemType::MARK_PREDECESSORS) ?  OFFSET_ELEMENTS : 0,
		};

		union {
			// Raking elements
			SizeT 							expand_raking_elements[RakingExpandGrid::TOTAL_RAKING_ELEMENTS];

			// Scratch elements
			struct {
				SizeT 						offset_scratch[OFFSET_ELEMENTS];
				VertexId 					predecessor_scratch[PARENT_ELEMENTS];
			};

			volatile VertexId 				warp_hashtable[WARPS][WARP_HASH_ELEMENTS];
			VertexId 						cta_hashtable[HASH_ELEMENTS];
		};
	};

	enum {
		THREAD_OCCUPANCY				= B40C_SM_THREADS(CUDA_ARCH) >> LOG_THREADS,
		SMEM_OCCUPANCY					= B40C_SMEM_BYTES(CUDA_ARCH) / sizeof(SmemStorage),
		CTA_OCCUPANCY  					= B40C_MIN(_MIN_CTA_OCCUPANCY, B40C_MIN(B40C_SM_CTAS(CUDA_ARCH), B40C_MIN(THREAD_OCCUPANCY, SMEM_OCCUPANCY))),

		VALID							= (CTA_OCCUPANCY > 0),
	};
};


} // namespace expand_contract_atomic
} // namespace bfs
} // namespace graph
} // namespace b40c

