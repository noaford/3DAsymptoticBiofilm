#ifndef _SYSTEMSOLVER_HPP_08172004_
#define _SYSTEMSOLVER_HPP_08172004_

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Linear system solver. 
/// \ingroup linalg
///
/// \todo Add documentation
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "NonConstructible.hpp"
#include "LinearAlgebra.hpp"

#include <stack>
#include <string>

namespace LinearAlgebra
{
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Static linear system solver class  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class DirectSolver : public NonConstructible
  {
  public:
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  /// \brief Performs a direct sparse matrix solve for \f$ A x = b \f$ 
	  /// \param A [in] Sparse matrix A
	  /// \param b [in] Sparse vector b
	  /// \param x [out] Sparse vector solution x
    /// \param Output [in] Set to true to write status messages. Default is true.    
	  /// \return Result of solve - true for success, false if solve failed.
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    static bool Solve( const LinearAlgebra::SparseMatrix& A, const LinearAlgebra::SparseVector& b, LinearAlgebra::SparseVector& x, bool Output = true );    
           
  private:
    typedef std::pair<int, int> SwapPair;
    typedef std::stack<SwapPair> Permutations;
        
    static void AssembleArrays( const LinearAlgebra::SparseMatrix& matA, const LinearAlgebra::SparseVector& vB, size_t& n, int*& ia, int*& ja, double*& A, double*& b, double*& x );
    static void AssembleVector( LinearAlgebra::SparseVector& vX, double* x );
        
    static void ShowDirectError( const int error );
        
    static std::string DirectSolverHeader;    
  }; 
}

#endif
