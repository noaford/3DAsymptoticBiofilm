#ifndef _LU_HPP_08302006_
#define _LU_HPP_08302006_

#include "CountedPtr.hpp"
#include "NonCopyable.hpp"
#include "SparseMatrix.hpp"
#include "Vector.hpp"


namespace LinearAlgebra
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief A %LU factorization of a sparse matrix.
  /// \ingroup linalg
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class LU : public NonCopyable
  {
  public:
  //@{
    /// \name Constructor and Destructor
    LU( const SparseMatrix& A );
    ~LU( void );
  //@}
    
  //@{
    /// \name Solver
    bool Solve( const Vector& b, Vector& x ) const;
  //@}

  private:
    void CreateMatrix( const SparseMatrix& A );
    void CreateVectors( const Vector& b, double*& pB, double*& pX ) const;
    void ShowError( const int error ) const;
  
    void* mPt[64];
    
    int mIparm[64];
    
    size_t mN;
    
    double* mA;
    int* mIa;
    int* mJa;
  };
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Counted pointer to a LU factorization.
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  typedef CountedPtr<LU> PLU;
}



#endif /*LU_HPP_*/
