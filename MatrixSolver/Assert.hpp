#ifndef _ASSERT_HPP_10262004_
#define _ASSERT_HPP_10262004_

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Expanded assert function.
/// \ingroup debug
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
#ifndef NASSERT
//------------------------------------------------------------------------------------------------------------------------------//
#include <string> 
#include <iostream>
#include <cstdlib>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Modified Assert Function.
///
/// Stops execution and displays an error message if the assertion fails. It is disabled if NASSERT if defined.
///
/// \param Assertion [in] Halts program if assertion is false, otherwise does nothing.
/// \param ErrorMsg [in] Error message to display if assertion is false.
/// \param File [in] Filename to display if assertion is false. Can be obtained from the preprocessor macro __FILE__.
/// \param Function [in] Function name to display if assertion is false. Can be obtained from the preprocessor macro __FUNCTION__.
/// \param Line [in] Line number to display if assertion is false. Can be obtained from the preprocessor macro __LINE__.
/// \ingroup debug
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline void Assert( const bool Assertion, const char* ErrorMsg, const char* File, const char* Function, const int Line )
{
  if ( !Assertion )
  { 
    std::cerr << "\n**** Assertion Failed **** " << ErrorMsg << " - File: " << File << " Function: " << Function << " Line: "; 
    std::cerr << Line << ".\n" << std::endl; 
    abort();
  }
}

#else
//------------------------------------------------------------------------------------------------------------------------------//
#define Assert( a, b, c, d, e )
//------------------------------------------------------------------------------------------------------------------------------//
#endif

#endif
