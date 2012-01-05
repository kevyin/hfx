
#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/set_operations.h>
#include "algorithms.h"
#include "utils.h"
#include "ion_series.h"

#include <stdint.h>

//#define DEBUG

/*
 * to lower
 */
__host__ __device__
uint8_t toLower(uint8_t up)
{
    return up + 32;
}

/*
 * to upper
 */
__host__ __device__
uint8_t toUpper(uint8_t lo)
{
    return lo - 32;
}



template <class Tin, class Tout>
struct residualByIdx : public thrust::unary_function<Tin, Tout>
{
    const Tout *d_residual_raw;
    const Tout delta;

    __host__ __device__
    residualByIdx(const Tout *_res, const Tout _d) : d_residual_raw(_res), delta(_d) {}

    __host__ __device__ Tout operator() (Tin idx)
    {
        return d_residual_raw[idx] + delta;
    }
};

// for a string of amino acids, apply modification specified by an unranked encoding
// return modified string of amino acids with modified acids(ma) in lower case
std::vector<uint8_t> 
genMIons
(
    std::vector<uint8_t>    ions, 
    std::vector<uint32_t>   unrank, 
    std::vector<uint8_t>    h_mod_ma, 
    std::vector<uint8_t>    h_mod_ma_count
)
{
    std::vector<uint8_t> mions; // modified ion str to return
    std::vector<uint32_t> cur_pep_ma_count(h_mod_ma.size());
    for (uint32_t i = 0; i < ions.size(); i++) {
        uint8_t acid = ions[i];
        // check if acid is a modable acid
        bool isMa = false;
        uint32_t ma_idx; // defined only if isMa == true
        for (uint32_t j = 0; j < h_mod_ma.size(); j++) {
            if (acid == h_mod_ma[j]) {
                isMa = true;
                ma_idx = j;    
                break;
            }
        }

        if (isMa) {
            // first record this finding and determine the ith 
            // occurrence of this ma
            uint32_t ith_ma = (cur_pep_ma_count[ma_idx]++);
            
            // check if this ma is to be modded according to unrank
            // first find the start to the relevant section of unrank
            uint32_t unrank_idx = 0;
            for (uint32_t k = 0; k < ma_idx; k++) {
                unrank_idx += h_mod_ma_count[k];                
            }
            bool modded = false;
            // check through unrank to see if the ith_ma is to be modded
            for (uint32_t l = 0; l < h_mod_ma_count[ma_idx]; l++) {
                if (ith_ma == unrank[unrank_idx + l]) {
                    mions.push_back(toLower(acid));
                    modded = true;
                    break;
                }
            }
            if (!modded)
                mions.push_back(acid);
        } else {
            mions.push_back(acid);
        }
    }
    return mions;
}

/**
 * Generate spectrums serially
 * also 
 */
void
getSpecNonParallel(
    uint32_t            *d_out_check_spec, 
    const uint8_t       *d_in_mions, 
    const float         *d_residual_raw, 
    const float         *d_mass_raw, 
    const uint8_t       *d_ions_raw, 
    const uint32_t      *d_tc_raw,
    const uint32_t      *d_tn_raw,
    const uint32_t      *d_mpep_idx_raw,
    const uint32_t      *d_mpep_unrank_raw,
    const uint32_t      num_mpep,
    const uint8_t       *d_mod_ma_raw,
    const uint8_t       *d_mod_ma_count_raw,
    const float         *d_mod_ma_mass_raw,
    const uint32_t      mod_num_ma,
    const float         res_delta,

    const uint32_t      max_charge,
    const uint32_t      len_spec
    
)
{
    // create a custom set of arrays to call addIons with
    // 
    // d_check_residual  mpep residual mass
    // d_check_mass      mass of amino acid
    // d_check_ions      amino acid chars. modified acid mass is lower case
    // d_check_tc        begin() of peptide
    // d_check_tn        end() of peptide
    // d_check_idx       "candidate" idx to above arrays. In this case should be every peptide
    // check_num_idx = num_mpep

    thrust::device_ptr<const float>         d_residual(d_residual_raw); 
    thrust::device_ptr<const float>         d_mass(d_mass_raw); 
    thrust::device_ptr<const uint8_t>       d_ions(d_ions_raw); 
    thrust::device_ptr<const uint32_t>      d_tc(d_tc_raw);
    thrust::device_ptr<const uint32_t>      d_tn(d_tn_raw);
    thrust::device_ptr<const uint32_t>      d_mpep_idx(d_mpep_idx_raw);
    thrust::device_ptr<const uint32_t>      d_mpep_unrank(d_mpep_unrank_raw);

    thrust::device_ptr<const uint8_t>       d_mod_ma(d_mod_ma_raw);
    thrust::device_ptr<const uint8_t>       d_mod_ma_count(d_mod_ma_count_raw);
    thrust::device_ptr<const float>         d_mod_ma_mass(d_mod_ma_mass_raw);
    

    // 
    std::vector<uint8_t> h_mions;
    std::vector<uint32_t> h_tc;
    std::vector<uint32_t> h_tn;

    uint32_t total_ma = thrust::reduce(d_mod_ma_count, d_mod_ma_count + mod_num_ma);
    std::vector<uint8_t> h_mod_ma(mod_num_ma);
    thrust::copy(d_mod_ma, d_mod_ma + mod_num_ma, h_mod_ma.begin());
    std::vector<uint8_t> h_mod_ma_count(mod_num_ma);
    thrust::copy(d_mod_ma_count, d_mod_ma_count + mod_num_ma, h_mod_ma_count.begin());

    // foreach unrank
    for (uint32_t i = 0; i < num_mpep; ++i) {
        // add mions and annotate tc tn
        uint32_t idx = d_mpep_idx[i];
        std::vector<uint8_t> ions(d_tn[idx] - d_tc[idx]);
        thrust::copy(d_ions + d_tc[idx], d_ions + d_tn[idx], ions.begin());

        std::vector<uint32_t> unrank(total_ma);
        thrust::copy(d_mpep_unrank + (i * total_ma), d_mpep_unrank + ((i+1) * total_ma), unrank.begin());

#ifdef _DEBUG
        for (std::vector<uint32_t>::iterator it = unrank.begin(); it != unrank.end(); it++) {
            std::cout << *it << " " ;
        }
        std::cout << " | ";
#endif
        std::vector<uint8_t> mions = genMIons(ions, unrank, h_mod_ma, h_mod_ma_count);

#ifdef _DEBUG
        for (std::vector<uint8_t>::iterator it = mions.begin(); it != mions.end(); it++) {
            std::cout << *it << " " ;
        }
        std::cout << std::endl;
#endif
        h_mions.insert(h_mions.end(), mions.begin(), mions.end());
        
        // add tc and tn's
        if (h_tc.size() == 0) {
            h_tc.push_back(0);
            h_tn.push_back(d_tn[idx] - d_tc[idx]);
        } else {
            h_tc.push_back(h_tn.back());
            h_tn.push_back(h_tn.back() + d_tn[idx] - d_tc[idx]);
        }

    }

    // fill out d_check_ arrays
    // d_check_residual  mpep residual mass
    thrust::device_vector<float> d_check_residual(num_mpep);
    thrust::transform(d_mpep_idx, d_mpep_idx + num_mpep, d_check_residual.begin(), residualByIdx<uint32_t, float>(d_residual_raw, res_delta));
    // d_check_mass      mass of amino acid
    
    thrust::device_vector<float>  d_check_mass(58);       // for ascii A to z
    thrust::copy(d_mass, d_mass + 26, d_check_mass.begin());
    for (uint32_t i = 0; i < mod_num_ma; ++i) {
        d_check_mass[toLower(d_mod_ma[i]) - 'A'] = d_mass[d_mod_ma[i] - 'A'] + d_mod_ma_mass[i];
    }
    // d_check_ions      amino acid chars. modified acid mass is lower case
    thrust::device_vector<uint8_t>  d_check_ions(h_mions.begin(), h_mions.end());
    // d_check_tc        begin() of peptide
    thrust::device_vector<uint32_t> d_check_tc(h_tc.begin(), h_tc.end());
    // d_check_tn        end() of peptide
    thrust::device_vector<uint32_t> d_check_tn(h_tn.begin(), h_tn.end());
    // d_check_idx       "candidate" idx to above arrays. In this case should be every peptide
    thrust::counting_iterator<uint32_t> first(0);
    thrust::counting_iterator<uint32_t> last = first + num_mpep;
    thrust::device_vector<uint32_t> d_check_idx(first, last);
    // check_num_idx = num_mpep
    uint32_t check_num_idx = num_mpep;


    ////////////////////////////////////////////////////////////////
#ifdef _DEBUG
    // DEBUG check the arrays
    // print some things
    std::cout << "mod_num_ma " << mod_num_ma << std::endl;
    std::cout << "res_delta " << res_delta << std::endl;
    std::cout << "max_charge " << max_charge << std::endl;
    std::cout << "len_spec " << len_spec << std::endl;
    
    // check residual with delta is correct
    std::cout << "checking residuals" << std::endl;
    for (uint32_t i = 0; i < num_mpep; ++i) {
        float p_mass = (d_residual[d_mpep_idx[i]] + res_delta);
        //std::cerr << "orig " << d_residual[d_mpep_idx[i]] << std::endl;
        //std::cerr << "p_mass " << p_mass << std::endl;
        //std::cerr << "s_mass " << d_check_residual[i] << std::endl;
        if (d_check_residual[i] != p_mass ){
            std::cerr << "orig " << d_residual_raw[d_mpep_idx[i]] << std::endl;
            std::cerr << d_check_residual[i] << " != " << p_mass << std::endl;
            exit(1);
        }
    }
    
    // d_mass
    for (uint8_t i = 'a'; i <= 'z'; ++i) {
        std::cout << "d_check_mass " << i << " " << d_check_mass[i - 'A'] << std::endl;
    }

    // d_check_ions      amino acid chars. modified acid mass is lower case
    thrust::device_ptr<const uint8_t> d_in_mions_th(d_in_mions);

    for (uint32_t i = 0; i < num_mpep; i++) {

        std::cout << "d_check_idx " << d_check_idx[i] << std::endl;
        uint32_t ch_c_idx = d_check_tc[i];
        uint32_t ch_n_idx = d_check_tn[i];

        std::cout << "check c " << ch_c_idx << " n " << ch_n_idx << std::endl;

        std::cout << "modded check  ";
        for (uint32_t j = ch_c_idx; j < ch_n_idx; ++j) {
            std::cout << d_check_ions[j] << " ";
        }
        std::cout << std::endl;

        std::cout << "from parallel ";
        for (uint32_t j = 0; j < ch_n_idx - ch_c_idx; ++j) {
            std::cout << d_in_mions_th[i*MAX_PEP_LEN + j] << " ";
        }
        std::cout << std::endl;
    }
#endif
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // call addions 
    addIons(d_out_check_spec,
            d_check_residual.data().get(),
            d_check_mass.data().get(),
            d_check_ions.data().get(),
            d_check_tc.data().get(),
            d_check_tn.data().get(),
            d_check_idx.data().get(),
            check_num_idx,
            max_charge,
            len_spec);
}
#undef DEBUG
