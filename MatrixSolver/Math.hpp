#ifndef _MATH_HPP_09122005_
#define _MATH_HPP_09122005_

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief %Math include file that allows the choice of math libraries.
///
/// Currently uses the Intel math libary when using the Intel C/C++ compiler and the regular C math libary for other compilers.
/// Not all functions are implemented when using the regular C math instead of the Intel.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "MathTypes.hpp"

#ifdef __INTEL_COMPILER
// Using Intel math library
  #include <mathimf.h>
#else
// Using standard C math libary
  #include <cmath>
  
  #include "Assert.hpp"
  
  // double cot( const double x )
  // { return( 1.0/tan(x) ); }
  // 
  // double sincos( const double x, double* sinval, double* cosval )
  // {
    // (*sinval) = sin(x);
    // (*cosval) = cos(x); 
  // }
  // 
  // double acosh( const double x )
  // {
    // Assert( false, "Not implemented", __FILE__, __FUNCTION__, __LINE__ );
    // return( 0.0 );
  // }
  // 
  // double asinh( const double x )
  // {
    // Assert( false, "Not implemented", __FILE__, __FUNCTION__, __LINE__ );
    // return( 0.0 );
  // }
  // 
  // double atanh( const double x )
  // {
    // Assert( false, "Not implemented", __FILE__, __FUNCTION__, __LINE__ );
    // return( 0.0 );
  // }
  // 
  // double sinhcosh( const double x, double* sinhval, double* coshval )
  // {
    // (*sinhval) = sinh(x);
    // (*coshval) = cosh(x);
  // }
  // 
  // double cbrt( const double x )
  // {
    // Assert( false, "Not implemented", __FILE__, __FUNCTION__, __LINE__ );
    // return( 0.0 );
  // }
  // 
  // double exp10( const double x )
  // { return( pow(10, x) ); }
  // 
  // double exp2( const double x )
  // { return( pow(2, x) ); }
  // 
  // double expm1( const double x )
  // { return( exp(x) - 1.0 ); }
  // 
  // double frexp( const double x, int* exp )
  // {
    // Assert( false, "Not implemented", __FILE__, __FUNCTION__, __LINE__ );
    // return( 0.0 );
  // }
  // 
  // double hypot( const double x, const double y )
  // { return( sqrt( x*x + y*y ) ); }
  // 
  // double invsqrt( const double x )
  // { return( 1.0/sqrt(x) ); }
// 
  // int ilogb( double x )
  // {
    // Assert( false, "Not implemented", __FILE__, __FUNCTION__, __LINE__ );
    // return( 0.0 );
  // }
  // 
  // double log1p( const double x )
  // { return( log( x + 1.0 ) ); }
  // 
  // double log2( const double x )
  // {
    // Assert( false, "Not implemented", __FILE__, __FUNCTION__, __LINE__ );
    // return( 0.0 );
  // }
  // 
  // double logb( double x )
  // {
    // Assert( false, "Not implemented", __FILE__, __FUNCTION__, __LINE__ );
    // return( 0.0 );
  // } 
 // 
  
#endif

#endif
