#ifndef _DIAGONALMATRIX_HPP_10122006_
#define _DIAGONALMATRIX_HPP_10122006_

#include "Assert.hpp"
#include "LinearAlgebra.hpp"

#include <ostream>

namespace LinearAlgebra
{
  class DiagonalMatrix
  {
  public:
    DiagonalMatrix( void );
    explicit DiagonalMatrix( const size_t N );    
    DiagonalMatrix( const size_t N, const double Value );
    DiagonalMatrix( const DiagonalMatrix& M );
    ~DiagonalMatrix( void );
    
    void SetValues( const double Value );
    void ZeroValues( void );
    
    DiagonalMatrix& operator= ( const DiagonalMatrix& M );
    
    DiagonalMatrix& operator+= ( const DiagonalMatrix& M );
    DiagonalMatrix& operator-= ( const DiagonalMatrix& M );
    DiagonalMatrix& operator*= ( const DiagonalMatrix& M );
    DiagonalMatrix& operator/= ( const DiagonalMatrix& M );
    
    DiagonalMatrix& operator*= ( const double c );
    DiagonalMatrix& operator/= ( const double c );

    size_t GetSize( void ) const;
    
    double operator() ( const int i ) const
    { 
      Assert( i < mSize, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      return( mData[i] ); 
    }
    
    double& operator() ( const int i )
    { 
      Assert( i < mSize, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      return( mData[i] ); 
    }
    
    friend std::ostream& operator<< ( std::ostream& o, const DiagonalMatrix& M );
    
    friend DiagonalMatrix operator+ ( const DiagonalMatrix& A, const DiagonalMatrix& B );
    friend DiagonalMatrix operator- ( const DiagonalMatrix& A, const DiagonalMatrix& B );
    friend DiagonalMatrix operator* ( const DiagonalMatrix& A, const DiagonalMatrix& B );
    friend DiagonalMatrix operator/ ( const DiagonalMatrix& A, const DiagonalMatrix& B );
    friend DiagonalMatrix operator- ( const DiagonalMatrix& A );
    
    friend DiagonalMatrix operator* ( const DiagonalMatrix& A, const double c );
    friend DiagonalMatrix operator/ ( const DiagonalMatrix& A, const double c );
    
    friend DiagonalMatrix operator* ( const double c, const DiagonalMatrix& B );
    
    friend Vector operator* ( const DiagonalMatrix& A, const Vector& x );
    friend Vector operator/ ( const DiagonalMatrix& A, const Vector& x );

    friend DiagonalMatrix inverse( const DiagonalMatrix& A );    

    friend double trace( const DiagonalMatrix& A );
    friend double frobnorm( const DiagonalMatrix& A );
    
    friend double absmax( const DiagonalMatrix& A );
    friend double absmin( const DiagonalMatrix& A );
    friend double max( const DiagonalMatrix& A );
    friend double min( const DiagonalMatrix& A );
    
  private:
    void AllocateMemory( void );
    void ReleaseMemory( void );
    void StraightCopy( const LinearAlgebra::DiagonalMatrix& M );
    void AllocateCopy( const LinearAlgebra::DiagonalMatrix& M );
    
    bool SizesMatch( const size_t N ) const;
    
    size_t mSize;
    double* mData;    
  };
}

#endif

