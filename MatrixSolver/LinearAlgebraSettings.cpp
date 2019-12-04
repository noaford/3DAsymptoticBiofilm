#include "LinearAlgebraSettings.hpp"

int LinearAlgebra::Settings::DirectSolverPerturb           = 13;
int LinearAlgebra::Settings::DirectSolverFileOutput        = 0;

void LinearAlgebra::Settings::GetSettings( ConfigFile& file )
{
  DirectSolverPerturb           = file.GetInteger( "DirectSolverPerturb", DirectSolverPerturb );
  DirectSolverFileOutput        = file.GetInteger( "DirectSolverFileOutput", DirectSolverFileOutput ); 
}

void LinearAlgebra::Settings::LoadSettings( const std::string Filename )
{  
  ConfigFile file( Filename );
  GetSettings( file );
}

void LinearAlgebra::Settings::InsertSettings( ConfigFile& file )
{
  file.SetInteger( "DirectSolverPerturb", DirectSolverPerturb );
  file.SetInteger( "DirectSolverFileOutput", DirectSolverFileOutput );  
}

void LinearAlgebra::Settings::SaveSettings( const std::string Filename )
{
  ConfigFile file( Filename );
  InsertSettings( file );
  file.SaveData();
}

void LinearAlgebra::Settings::DisplaySettings( void )
{
  std::cout << "Linear Algebra Settings:" << std::endl;
  std::cout << "----------------------------------------------------------" << std::endl;
  std::cout << "  DirectSolver -->  Amount to Perturb Pivots = 1e-" << DirectSolverPerturb << std::endl;  
  std::cout << "  DirectSolver --> Output Solve Data to File = " << DirectSolverFileOutput << std::endl;   
  std::cout << "----------------------------------------------------------" << std::endl;
}

