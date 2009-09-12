/*
 * Module    : IonSeries
 * Copyright : (c) 2009 Trevor L. McDonell
 * License   : BSD
 */

#include "mass.h"
#include "utils.h"
#include "kernels.h"


/*
 * Convert a given mass into a mass/charge ratio, and locate the appropriate
 * spectrum bin for the peak.
 */
__device__ int
binIonMZ(float mass, int charge)
{
    int bin = rintf((mass + massH * charge) / (charge * binWidthMono));
#ifdef __DEVICE_EMULATION__
    assert(bin >= 0 && bin < 2048);
#endif
    return bin;
}


/*
 * Add a spectral peak for each fragment ion location, as well as the peaks
 * corresponding to the neutral losses of H2O and NH3.
 */
__device__ void
addIonsAB(float mass, int charge, int *spec)
{
    // A
    atomicMax(&spec[binIonMZ(mass - massCO, charge)], 10);

    // B
    int m = binIonMZ(mass, charge);

    atomicMax(&spec[m],   50);
    atomicMax(&spec[m+1], 25);
    atomicMax(&spec[m-1], 25);
    atomicMax(&spec[binIonMZ(mass - massH2O, charge)], 10);
    atomicMax(&spec[binIonMZ(mass - massNH3, charge)], 10);
}


__device__ void
addIonsY(float mass, int charge, int *spec)
{
    int m = binIonMZ(mass + massH2O, charge);

    atomicMax(&spec[m],   50);
    atomicMax(&spec[m+1], 25);
    atomicMax(&spec[m-1], 25);
    atomicMax(&spec[binIonMZ(mass - massNH3, charge)], 10);
}


/*
 * Add a spectral peak for each fragment ion location. The output spectrum array
 * must exist and be initialised to zero.
 */
template <bool lengthIsPow2>
__global__ static void
addIons_core
(
    int          max_charge,
    float        *b_ions,
    float        *y_ions,
    int          *spec,
    unsigned int len_ions
)
{
    unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;

    if (lengthIsPow2 || idx < len_ions)
    {
        int   charge = 1;
        float b_mass = b_ions[idx];
        float y_mass = y_ions[idx];

        do
        {
            addIonsAB(b_mass, charge, spec);
            addIonsY (y_mass, charge, spec);
        }
        while (++charge < max_charge);
    }
}


void
addIons
(
    int          max_charge,
    float        *b_ions,
    float        *y_ions,
    int          *spec,
    unsigned int len_ions,
    unsigned int len_spec
)
{
    unsigned int threads = min(ceilPow2(len_ions), 512);
    unsigned int blocks  = (len_ions + threads - 1) / threads;

    (void) len_spec;

    if (isPow2(len_ions))
        addIons_core<true><<<blocks,threads>>>(max_charge, b_ions, y_ions, spec, len_ions);
    else
        addIons_core<false><<<blocks,threads>>>(max_charge, b_ions, y_ions, spec, len_ions);
}

