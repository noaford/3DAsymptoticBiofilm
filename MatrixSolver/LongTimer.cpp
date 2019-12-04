#include "LongTimer.hpp"

#include <ostream>
#include <sstream>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Generates a string in Hrs:Min:Sec format using the time elapsed.
/// \return String representing the time elapsed in HH:MM:SS formate.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string LongTimer::HrsMinSec( void ) const
{
    time_t sec = Elapsed();
    
    unsigned long hours = static_cast<unsigned long>( sec/( 60*60) );
    unsigned long min   = static_cast<unsigned long>( sec/60 ) - hours*60;
    sec -= ( hours*60*60  + min*60 );
    
    std::ostringstream oss;
    
    oss << hours << ":";
    
    if (min < 10)     // Insert a leading zero if minutes is a single digit 
      oss << "0";
    
    oss << min << ":";
    
    if (sec < 10)     // Insert a leading zero if seconds is a single digit
      oss << "0";
    
    oss << sec;
    
    return( oss.str() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Generates a message string from a given message and the timer's value
/// \return String generated from the given message and the timer's value.  It is of the form "'message' in HH::MM::SS seconds."
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string LongTimer::Message( const std::string message ) const
{
  std::ostringstream oss;
  oss << message << " in " << HrsMinSec() << " seconds.";
  
  return( oss.str() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Stream output operator.
/// \return Output stream with the time elapsed inserted in Hrs:Min:Sec format.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::ostream& operator<< ( std::ostream& o, const LongTimer& t )
{    
  return( o << t.HrsMinSec() );
}
