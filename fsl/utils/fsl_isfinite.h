#ifndef fsl_isfinite_h
#define fsl_isfinite_h

#include <math.h>

/* Cross-platform isfinite support */
#if defined(sun) || defined(__sun)
#include <ieeefp.h>
#define isfinite(_fsl_n) (finite(_fsl_n) && !isnan(_fsl_n))
#else
#include <float.h>
#define isfinite(_fsl_n) (_finite(_fsl_n) && !_isnan(_fsl_n))
#endif

#endif

