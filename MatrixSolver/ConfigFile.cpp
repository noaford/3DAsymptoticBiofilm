#include "ConfigFile.hpp"

std::string ConfigFile::mError[4] = { "Missing '='", "Blank key", "Blank value", "Duplicate key" };
  
void ConfigFile::LoadData( void )
{
  std::ifstream file( mFilename.c_str() );
  
  if ( !file )
  { // File not opened.
    std::cerr << "Cannot open configuration file '" << mFilename << "'." << std::endl;
    return;
  }
  
  std::string strOriginal, strLine;
  int nLine = 1;
  
  while ( getline(file, strOriginal) )       // Get a new line from the file.
  { 
    strLine = strOriginal;
    StripComment( strLine );
    
    if ( !IsBlank( strLine ) )
    {
      std::string::size_type idx = strLine.find("=");
      if ( idx == std::string::npos )
        ShowError( nLine, mError[0], mFilename, strOriginal );
      
      std::string Key, Value;      
      TokenizeLine( strLine, idx, Key, Value );
      
      if ( IsBlank( Key ) )
        ShowError( nLine, mError[1], mFilename, strOriginal );
      
      if ( IsBlank( Value ) )
        ShowError( nLine, mError[2], mFilename, strOriginal );
      
      if ( mKeyMap.find(Key) != mKeyMap.end() )
        ShowError( nLine, mError[3], mFilename, strOriginal );
      else
        mKeyMap.insert( KeyPair(Key, Value) );
    }
    ++nLine;
  }
}

