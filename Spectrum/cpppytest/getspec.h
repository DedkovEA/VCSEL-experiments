/*
   Declares the functions to be imported by our application, and exported by our
   DLL, in a flexible and elegant way.
*/

#include <complex>

/* You should define ADD_EXPORTS *only* when building the DLL. */
#ifdef ADD_EXPORTS
    #define ADDAPI __declspec(dllexport)
#else
    #define ADDAPI __declspec(dllimport)
#endif

/* Define real type for all calculations. */
#ifndef INTRINSIC_TYPE
    #define INTRINSIC_TYPE double
#endif
typedef INTRINSIC_TYPE floating;
typedef std::complex<floating> cfloating;

/* Define calling convention in one place, for convenience. */
#define ADDCALL __cdecl

/* Make sure functions are exported with C linkage under C++ compilers. */

#ifdef __cplusplus
extern "C"
{
#endif

/* Declare functions using the above definitions. */


// Some memory managment functions for python

ADDAPI cfloating* ADDCALL reserve_complex_array(int size);
ADDAPI floating* ADDCALL reserve_array(int size);

ADDAPI cfloating* ADDCALL reserve_twists(int Npow);
ADDAPI unsigned int* ADDCALL reserve_bitrev(int Npow);


ADDAPI void ADDCALL free_complex_array(cfloating* ptr);
ADDAPI void ADDCALL free_array(floating* ptr);
ADDAPI void ADDCALL free_uint_array(unsigned int* ptr);


ADDAPI void ADDCALL sdeeval(floating* specx, floating* specy, cfloating* Ex, cfloating* Ey, cfloating* tmpEx, cfloating* tmpEy,  // arrays for SDE samples and output spectras 
                            cfloating* twist, unsigned int* bitrev, int Npow, int skip, int Nav, int tauDt, floating Dt,            // some necessary variables
                            floating alpha, floating kappa, floating gamma, floating gamma_d, floating gamma_a,                  // parameters
                            floating gamma_p, floating beta, floating mu, floating C_sp, floating N_th, floating N_tr);


#ifdef __cplusplus
} // __cplusplus defined.
#endif