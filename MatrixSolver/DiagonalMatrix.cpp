#include "DiagonalMatrix.hpp"

#include "Assert.hpp"

#include "Vector.hpp"

LinearAlgebra::DiagonalMatrix::DiagonalMatrix( void ) : mSize(1), mData(NULL) 
{ AllocateMemory(); }

LinearAlgebra::DiagonalMatrix::DiagonalMatrix( const size_t N ) : mSize(N), mData(NULL) 
{ AllocateMemory(); }

LinearAlgebra::DiagonalMatrix::DiagonalMatrix( const size_t N, const double Value ) : mSize(N), mData(NULL)
{
  AllocateMemory();
  SetValues( Value );
}

LinearAlgebra::DiagonalMatrix::DiagonalMatrix( const LinearAlgebra::DiagonalMatrix& M ) : mSize(M.mSize), mData(NULL)
{
  AllocateMemory();
  StraightCopy( M ); 
}

LinearAlgebra::DiagonalMatrix::~DiagonalMatrix( void ) 
{ ReleaseMemory(); }

void LinearAlgebra::DiagonalMatrix::SetValues( const double Value )
{ std::fill( mData, mData+mSize, Value ); }

void LinearAlgebra::DiagonalMatrix::ZeroValues( void )
{ SetValues( 0.0 ); }

LinearAlgebra::DiagonalMatrix& LinearAlgebra::DiagonalMatrix::operator= ( const LinearAlgebra::DiagonalMatrix& M )
{
  if ( SizesMatch( M.mSize ) )
    StraightCopy( M );
  else
  {
    mSize = M.mSize,
    AllocateCopy( M );
  }
  
  return( *this );
}

LinearAlgebra::DiagonalMatrix& LinearAlgebra::DiagonalMatrix::operator+= ( const LinearAlgebra::DiagonalMatrix& M )
{
  Assert( SizesMatch( M.mSize ), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  std::transform( mData, mData+mSize, M.mData, mData, std::plus<double>() ); 
  return(*this);  
}

LinearAlgebra::DiagonalMatrix& LinearAlgebra::DiagonalMatrix::operator-= ( const LinearAlgebra::DiagonalMatrix& M )
{
  Assert( SizesMatch( M.mSize ), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  std::transform( mData, mData+mSize, M.mData, mData, std::minus<double>() ); 
  return(*this);  
}

LinearAlgebra::DiagonalMatrix& LinearAlgebra::DiagonalMatrix::operator*= ( const LinearAlgebra::DiagonalMatrix& M )
{
  Assert( SizesMatch( M.mSize ), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  std::transform( mData, mData+mSize, M.mData, mData, std::multiplies<double>() ); 
  return(*this);  
}

LinearAlgebra::DiagonalMatrix& LinearAlgebra::DiagonalMatrix::operator/= ( const LinearAlgebra::DiagonalMatrix& M )
{
  Assert( SizesMatch( M.mSize ), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  std::transform( mData, mData+mSize, M.mData, mData, std::divides<double>() ); 
  return(*this);  
}

LinearAlgebra::DiagonalMatrix& LinearAlgebra::DiagonalMatrix::operator*= ( const double c )
{
  std::transform( mData, mData+mSize, mData, std::bind2nd(std::multiplies<double>(),c) ); 
  return(*this);  
}

LinearAlgebra::DiagonalMatrix& LinearAlgebra::DiagonalMatrix::operator/= ( const double c )
{
  std::transform( mData, mData+mSize, mData, std::bind2nd(std::divides<double>(),c) ); 
  return(*this);  
}

size_t LinearAlgebra::DiagonalMatrix::GetSize( void ) const
{ return( mSize ); }

std::ostream& LinearAlgebra::operator<< ( std::ostream& o, const LinearAlgebra::DiagonalMatrix& M )
{
  for ( int i = 0; i < M.mSize; ++i )
    o << M.mData[i] << "\n";
    
  return( o );
}

LinearAlgebra::DiagonalMatrix LinearAlgebra::operator+ ( const LinearAlgebra::DiagonalMatrix& A, const LinearAlgebra::DiagonalMatrix& B )
{
  Assert( A.SizesMatch(B.mSize), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );  
  LinearAlgebra::DiagonalMatrix result( A );
  
  return( result += B );
}

LinearAlgebra::DiagonalMatrix LinearAlgebra::operator- ( const LinearAlgebra::DiagonalMatrix& A, const LinearAlgebra::DiagonalMatrix& B )
{
  Assert( A.SizesMatch(B.mSize), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  LinearAlgebra::DiagonalMatrix result( A );
  
  return( result -= B );
}

LinearAlgebra::DiagonalMatrix LinearAlgebra::operator* ( const LinearAlgebra::DiagonalMatrix& A, const LinearAlgebra::DiagonalMatrix& B )
{
  Assert( A.SizesMatch(B.mSize), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  LinearAlgebra::DiagonalMatrix result( A );
  
  return( result *= B );
}

LinearAlgebra::DiagonalMatrix LinearAlgebra::operator/ ( const LinearAlgebra::DiagonalMatrix& A, const LinearAlgebra::DiagonalMatrix& B )
{
  Assert( A.SizesMatch(B.mSize), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  LinearAlgebra::DiagonalMatrix result( A );
  
  return( result /= B );
}

LinearAlgebra::DiagonalMatrix LinearAlgebra::operator* ( const LinearAlgebra::DiagonalMatrix& A, const double c )
{
  LinearAlgebra::DiagonalMatrix result( A );
  
  return( result *= c );
}

LinearAlgebra::DiagonalMatrix LinearAlgebra::operator/ ( const LinearAlgebra::DiagonalMatrix& A, const double c )
{
  LinearAlgebra::DiagonalMatrix result( A );
  
  return( result /= c );
}
LinearAlgebra::DiagonalMatrix LinearAlgebra::operator* ( const double c, const LinearAlgebra::DiagonalMatrix& B )
{
  LinearAlgebra::DiagonalMatrix result( B );
  
  return( result *= c );
}

LinearAlgebra::Vector LinearAlgebra::operator/ ( const LinearAlgebra::DiagonalMatrix& A, const LinearAlgebra::Vector& x )
{
  const size_t nSize = x.GetSize();
  Assert( A.SizesMatch(nSize), "Vector and Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  LinearAlgebra::Vector result( x );

  for ( int i = 0; i < nSize; ++i )
    result(i) /= A(i);
  
  return( result );
}

LinearAlgebra::Vector LinearAlgebra::operator* ( const LinearAlgebra::DiagonalMatrix& A, const LinearAlgebra::Vector& x )
{
  const size_t nSize = x.GetSize();
  Assert( A.SizesMatch(nSize), "Vector and Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  LinearAlgebra::Vector result( x );

  for ( int i = 0; i < nSize; ++i )
    result(i) *= A(i);
  
  return( result );
}

LinearAlgebra::DiagonalMatrix LinearAlgebra::inverse( const LinearAlgebra::DiagonalMatrix& A )
{
  LinearAlgebra::DiagonalMatrix result( A.mSize );
  
  for ( int i = 0; i < A.mSize; ++i )
    result(i) = 1.0/A(i);
    
  return( result );
}

double LinearAlgebra::trace( const LinearAlgebra::DiagonalMatrix& A )
  { return( std::accumulate( A.mData, A.mData+A.mSize, 1.0, std::multiplies<double>() ) ); }

double LinearAlgebra::frobnorm( const LinearAlgebra::DiagonalMatrix& A )
{ return( std::inner_product( A.mData, A.mData+A.mSize, A.mData, 0.0 ) ); }

void LinearAlgebra::DiagonalMatrix::AllocateMemory( void )
{
  Assert( mSize > 0, "Negative Matrix Size", __FILE__, __FUNCTION__, __LINE__ );
  
  mData = new double[mSize];
  Assert( mData != NULL, "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );  
}

void LinearAlgebra::DiagonalMatrix::ReleaseMemory( void )
{ SafeDeleteArray( mData ); }

void LinearAlgebra::DiagonalMatrix::StraightCopy( const LinearAlgebra::DiagonalMatrix& M )
{ memcpy( mData, M.mData, sizeof(double)*mSize ); }

void LinearAlgebra::DiagonalMatrix::AllocateCopy( const LinearAlgebra::DiagonalMatrix& M )
{
  ReleaseMemory();
  AllocateMemory();
  StraightCopy( M );
}

bool LinearAlgebra::DiagonalMatrix::SizesMatch( const size_t nSize ) const
{ return( mSize == nSize ); }

LinearAlgebra::DiagonalMatrix LinearAlgebra::operator- ( const LinearAlgebra::DiagonalMatrix& A )
{ return( -1.0 * A ); }

// absmax was here
double LinearAlgebra::absmax( const LinearAlgebra::DiagonalMatrix& A )
{ 
  double Max = 0.0;
  for ( int i = 0; i < A.mSize; ++i )
  {
    const double Value =  fabs(A.mData[i]);
    Max = std::max( Max, Value );
  }
  
  return( Max ); 
}

double LinearAlgebra::absmin( const LinearAlgebra::DiagonalMatrix& A )
{
  double Min = std::numeric_limits<double>::max();
  for ( int i = 0; i < A.mSize; ++i )
  {
    const double Value =  fabs(A.mData[i]);
    Min = std::min( Min, Value );
  }
  
  return( Min );
}

double LinearAlgebra::max( const LinearAlgebra::DiagonalMatrix& A )
{ return( *std::max_element( A.mData, A.mData+A.mSize ) ); }

double LinearAlgebra::min( const LinearAlgebra::DiagonalMatrix& A )
{ return( *std::min_element( A.mData, A.mData+A.mSize ) ); }
