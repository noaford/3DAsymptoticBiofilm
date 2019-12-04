#ifndef _NONCONSTRUCTIBLE_HPP_10262004_
#define _NONCONSTRUCTIBLE_HPP_10262004_

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Base class to prevent class construction.
/// \ingroup base
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Inherit from this class to prevent the child class from being constructed.
/// \ingroup base
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class NonConstructible
{
protected:
  NonConstructible( void );  
};

#endif
