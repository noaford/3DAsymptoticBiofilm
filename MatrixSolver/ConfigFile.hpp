#ifndef _CONFIGFILE_HPP_03012005_
#define _CONFIGFILE_HPP_03012005_

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Configuration file manager.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "NonCopyable.hpp"

#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Class for loading and saving configuration information in a file.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class ConfigFile : public NonCopyable
{
public:
//@{
  /// \name Constructor
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Constructs a ConfigFile and loads the data from the given file.
  ///
  /// Creates a new file the given filename does not exist.
  /// \param Filename [in] Filename of the configuration file to load.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ConfigFile( const std::string& Filename ) : mFilename( Filename ) { LoadData(); }
  
  // Uses compiler generated destructor.
  // Uses compiler generated copy constructor.
  // Uses compiler generated assignment operator.
//@}

//@{
  /// \name File output functions
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Saves the data to the configuration file. 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SaveData( void ) const
  {
    std::ofstream file( mFilename.c_str() );
    for ( ConstKeyMapIter iter = mKeyMap.begin(); iter != mKeyMap.end(); ++iter )
      file << iter->first << " = " << iter->second << std::endl;   
  }
//@}

//@{
  /// \name Memory functions
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Releases the internal memory used for storing the settings. 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void ReleaseData( void ) { KeyMap().swap( mKeyMap ); }
//@}
  
//@{
  /// \name Data input functions.  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Gets an integer from the config file.
  /// \param Key [in] Key name in the config file.
  /// \param Default [in] Default value to return if the key is not found.
  /// \return Integer associated with the given key or default if the key is not found.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int GetInteger( const std::string& Key, const int Default = 0 ) const
  { 
    ConstKeyMapIter iter = mKeyMap.find( Key );
    return( ( iter != mKeyMap.end() ) ? atoi( iter->second.c_str() ) : Default );
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Gets a double from the config file.
  /// \param Key [in] Key name in the config file.
  /// \param Default [in] Default value to return if the key is not found.
  /// \return Double associated with the given key or default if the key is not found.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  double GetDouble( const std::string& Key, const double Default = 0.0 ) const
  {
    ConstKeyMapIter iter = mKeyMap.find( Key );
    return( ( iter != mKeyMap.end() ) ? atof( iter->second.c_str() ) : Default );
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Gets a boolean from the config file.
  /// \param Key [in] Key name in the config file.
  /// \param Default [in] Default value to return if the key is not found.
  /// \return Boolean associated with the given key or default if the key is not found.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool GetBoolean( const std::string& Key, const bool Default = false ) const
  {
    ConstKeyMapIter iter = mKeyMap.find( Key );
    return( ( iter != mKeyMap.end() ) ? atoi( iter->second.c_str() ) : Default );
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Gets a string from the config file.
  /// \param Key [in] Key name in the config file.
  /// \param Default [in] Default value to return if the key is not found.
  /// \return String associated with the given key or default if the key is not found.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::string GetString( const std::string& Key, const std::string& Default = std::string("") ) const
  {
    ConstKeyMapIter iter = mKeyMap.find( Key );
    return( ( iter != mKeyMap.end() ) ? iter->second : Default );
  }
//@}

//@{
  /// \name Data output functions
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Sets the value of an integer field in the config file.
  /// \param Key [in] Key name in the config file. If key does not already exist, then the key is added to the file.
  /// \param Value [in] Value to assign to the key.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SetInteger( const std::string& Key, const int Value )
  { 
    std::ostringstream ss;
    ss << Value;
    mKeyMap[Key] = ss.str();
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Sets the value of a double field in the config file.
  /// \param Key [in] Key name in the config file. If key does not already exist, then the key is added to the file.
  /// \param Value [in] Value to assign to the key.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SetDouble( const std::string& Key, const double Value )
  { 
    std::ostringstream ss;
    ss << Value;
    mKeyMap[Key] = ss.str();
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Sets the value of a boolean field in the config file.
  /// \param Key [in] Key name in the config file. If key does not already exist, then the key is added to the file.
  /// \param Value [in] Value to assign to the key.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SetBoolean( const std::string& Key, const bool Value )
  { 
    std::ostringstream ss;
    ss << Value;
    mKeyMap[Key] = ss.str();
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Sets the value of a string field in the config file.
  /// \param Key [in] Key name in the config file. If key does not already exist, then the key is added to the file.
  /// \param Value [in] Value to assign to the key.
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void SetString( const std::string& Key, const std::string& Value ) { mKeyMap[Key] = Value; }
//@}
  
private:
  void LoadData( void );
  
  void ShowError( const int Line, const std::string& Error, const std::string& Filename, const std::string& Code ) const
  {
    std::cerr << Error << " on line " << Line << " in " << Filename << std::endl;
    std::cerr << "\t" << Code << std::endl;
  }
  
  void StripComment( std::string& strLine ) const
  {
    std::string::size_type idx = strLine.find("#");   // Find the comment token (#)
    if ( idx != std::string::npos )
      strLine.erase(idx);
  }

  void TokenizeLine( const std::string& Line, const std::string::size_type idx, std::string& Key, std::string& Value ) const
  {  
    Key = Line.substr(0, idx-1);
    Value = Line.substr(idx+1);
  }
  
  void StripWhitespace( std::string& str ) const
  {
    std::string::size_type idx;
    while ( (idx = str.find(" ") ) != std::string::npos )
      str.erase( idx, 1 );
    
    while ( (idx = str.find("\t") ) != std::string::npos )
      str.erase( idx, 1 );
  }
  
  bool IsBlank( const std::string& str ) const
  {
    std::string strTemp( str );
    StripWhitespace( strTemp );
    return( strTemp.empty() );
  }

  
  typedef std::map< std::string, std::string > KeyMap;
  typedef KeyMap::value_type KeyPair;
  typedef KeyMap::iterator KeyMapIter;
  typedef KeyMap::const_iterator ConstKeyMapIter;
  
  std::string mFilename;
  KeyMap mKeyMap;
  
  static std::string mError[4];
};

#endif
