#ifndef _SAFEDELETE_HPP_10262004_
#define _SAFEDELETE_HPP_10262004_

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Safe memory deletion functions.
/// \ingroup mem
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <cstdlib>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Safely deletes a pointer.
/// \param p [in, out] Pointer to delete.  The pointer is set to NULL.
/// \ingroup mem
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class TYPE>
inline void SafeDelete( TYPE& p ) 
{ 
  if ( p != NULL ) 
  { 
    delete p; 
    p = NULL; 
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Safely deletes an array pointer.
/// \param p [in, out] Pointer to an array.  The pointer is set to NULL.
/// \ingroup mem
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class TYPE>
inline void SafeDeleteArray( TYPE& p )
{
  if ( p != NULL )
  {
    delete[] p;
    p = NULL;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Safely deletes an array of pointers.
/// \warning Does not delete the array memory, only the memory pointed to by the the array elements.  Must call SafeDeleteArray to
///          to delete the array afterwards.
/// \param p [in, out] Pointer to an array of pointers.  Each element of the array is set to NULL.
/// \param n [in] Number of pointers in the array.
/// \ingroup mem
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class TYPE>
inline void SafeDeletePointerArray( TYPE& p, const int n )
{
  if ( p != NULL )
  {
    for ( int i = 0; i < n; ++i ) 
      SafeDelete( p[i] );
    p = NULL;
  }
}

#endif
