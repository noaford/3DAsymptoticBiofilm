#include "Matrix.hpp"

#include "SparseMatrix.hpp"

#ifdef _USE_INTEL_VML_
  #include "mkl_vml.h"
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sparse matrix copy constructor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Matrix::Matrix( const SparseMatrix& m ) : mM(m.mN), mN(m.mN), mSize(m.mN*m.mN), mData(NULL) 
{ 
  AllocateMemory(); 
  StraightCopy(m); 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sparse matrix assignment operator.
///
/// If the source and the destination matrix sizes match, then the value are copied, otherwise memory is reallocated in the 
/// destination matrix to hold the values of the source matrix.
/// \param m [in] Source matrix.
/// \return Reference to the destination matrix after the assignment.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Matrix& LinearAlgebra::Matrix::operator= ( const SparseMatrix& m )
{ 
  if ( SizesMatch(m.mN, m.mN) ) 
    StraightCopy(m); 
  else 
  {
    mM = m.mN;
    mN = m.mN;
    AllocateCopy(m);
  }
  return( *this ); 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sparse matrix inplace addition operator.
///
/// \warning The size of the source and destination matrices must match.
/// \param m [in] Source matrix.
/// \return Reference to the destination matrix after the addition.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Matrix&	LinearAlgebra::Matrix::operator+= ( const LinearAlgebra::SparseMatrix& m )
{
  Assert( m.mN == mM, "Matrix Sizes do not Match", __FILE__, __FUNCTION__, __LINE__ );
  Assert( m.mN == mN, "Matrix Sizes do not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  for ( int i = 0; i < m.mN; ++i )
  {
    SparseMatrix::RowIterator End = m.mRows[i].end();
    for ( SparseMatrix::RowIterator iter = m.mRows[i].begin(); iter != End; ++iter )
      mData[ Index(i, iter->first) ] += iter->second;
  }
  
  return( *this );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sparse matrix inplace subtraction operator.
///
/// \warning The size of the source and destination matrices must match.
/// \param m [in] Source matrix.
/// \return Reference to the destination matrix after the subtraction.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Matrix&	LinearAlgebra::Matrix::operator-= ( const LinearAlgebra::SparseMatrix& m )
{  
  Assert( m.mN == mM, "Matrix Sizes do not Match", __FILE__, __FUNCTION__, __LINE__ );
  Assert( m.mN == mN, "Matrix Sizes do not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  for ( int i = 0; i < m.mN; ++i )
  {
    SparseMatrix::RowIterator End = m.mRows[i].end();
    for ( SparseMatrix::RowIterator iter = m.mRows[i].begin(); iter != End; ++iter )
      mData[ Index(i, iter->first) ] -= iter->second;
  }
  
  return( *this );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Dense matrix inplace multiplication operator.
///
/// \note The number of columns in the destination matrix will change if the source matrix is not square.
/// \warning The number of rows of the source matrix must equal the number of columns in the destination matrix.
/// \todo Optimize code for when the source matrix is square.
/// \param m [in] Source matrix.
/// \return Reference to teh destination matrix after the multiplication.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Matrix& LinearAlgebra::Matrix::operator*= ( const LinearAlgebra::Matrix& m )
{
  Assert( (mN == m.mM), "Incompatible Matrix Sizes", __FILE__, __FUNCTION__, __LINE__  );  
  double* pData = new double[mM*m.mN];
  
  memset( pData, 0, sizeof(double)*mSize );
  
  for ( int i = 0; i < mM; ++i )
    for ( int j = 0; j < m.mN; ++j )
	    for ( int k = 0; k < mN; ++k )
	      pData[ Index(i,j) ] += mData[ Index(i,k) ] * m.mData[ m.Index(k,j) ];

  std::swap(mData, pData);
  SafeDeleteArray(pData);
  
  mN = m.mN;
  
  return( *this );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sparse matrix inplace multiplication operator.
///
/// \note The number of columns in the destination matrix will change if the source matrix is not square.
/// \warning The number of rows of the source matrix must equal the number of columns in the destination matrix.
/// \todo Optimize code for when the source matrix is square.
/// \param m [in] Source matrix.
/// \return Reference to teh destination matrix after the multiplication.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Matrix& LinearAlgebra::Matrix::operator*= ( const SparseMatrix& m )
{
  Assert( SizesMatch(m.mN, m.mN), "Incompatible Matrix Sizes", __FILE__, __FUNCTION__, __LINE__  );
  double* pData = new double[mSize];
  Assert( pData, "Unable to allocate memory", __FILE__, __FUNCTION__, __LINE__  );
  memset( pData, 0, sizeof(double)*mSize );
  
  for ( int i = 0; i < mN; ++i )
    for ( int k = 0; k < mN; ++k )
    {
      SparseMatrix::RowIterator End = m.mRows[k].end();
      for ( SparseMatrix::RowIterator iter = m.mRows[k].begin(); iter != End; ++iter )
        pData[ Index(i,iter->first) ] += mData[ Index(i,k) ] * iter->second;
    }
  
  std::swap(mData, pData);
  SafeDeleteArray(pData);
    
  mM = m.mN;
  mN = m.mN;
  return( *this );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Assigns a matrix row the values of a given vector.
/// \warning The length of the vector must equal the number of columns in the matrix.
/// \param i [in] Row index.
/// \param v [in] Source vector.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::Matrix::SetRow( const int i, const Vector& v )
{
  Assert( v.mSize == mN, "Matrix/Vector Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  std::copy( v.mData, v.mData+v.mSize, mData+Index(i,0) ); 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Assigns a matrix column the values of a given vector.
/// \warning The length of the vector must equal the number of rows in the matrix.
/// \param j [in] Column index.
/// \param v [in] Source vector.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::Matrix::SetColumn( const int j, const Vector& v ) 
{     
  for ( int i = 0; i < mM; ++i ) 
    mData[Index(i,j)] = v(i); 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Inserts a sparse matrix.
///
/// Assign the values of the matrix to be
/// \f[ A_{m+i,n+j} = B_{m, n} \text{for} m, n \leq N\f]
/// where \f$ A \f$ is the destination matrix, \f$ B \f$ is the source matrix, and \f$ N \f$ is the size of the source matrix.
/// \warning \f$ N + i \f$ must be less than the number of rows in the matrix.
/// \warning \f$ N + j \f$ must be less than the number of columns in the matrix.
/// \param m [in] Source matrix.
/// \param i [in] Row offset.
/// \param j [in] Column offset.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::Matrix::InsertBlock( const LinearAlgebra::SparseMatrix& m, const int i, const int j )
{   
  for ( int I = 0; I < m.mN; ++I )
  {
    SparseMatrix::RowIterator End = m.mRows[i].end();
    for ( SparseMatrix::RowIterator iter = m.mRows[i].begin(); iter != End; ++iter )
      mData[ Index(i+I, j+(iter->first)) ] = (iter->second);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Inserts a dense matrix.
///
/// Assign the values of the matrix to be
/// \f[ A_{m+i,n+j} = B_{m, n} \text{for} m \leq M, n \leq N\f]
/// where \f$ A \f$ is the destination matrix, \f$ B \f$ is the source matrix, \f$ M \f$ is the number of rows and \f$ N \f$ is 
/// the number of columns in the source matrix.
/// \warning \f$ M + i \f$ must be less than the number of rows in the matrix.
/// \warning \f$ N + j \f$ must be less than the number of columns in the matrix.
/// \param m [in] Source matrix.
/// \param i [in] Row offset.
/// \param j [in] Column offset.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::Matrix::InsertBlock( const LinearAlgebra::Matrix& m, const int i, const int j )
{
//  std::cout << "Insert Block:\n"; 
//  std::cout << "          Matrix Size = " << mM << "x" << mN << "\n";
//  std::cout << "      Insertion Index = " << i << "x" << j << "\n";
//  std::cout << " Inserted Matrix Size = " << m.mM << "x" << m.mN << "\n";
//  std::cout << "  Insertion End Index = " << m.mM + i << "x" << m.mN + j << "\n";
  
  Assert( mM <= (m.mM+i), "Inserted Matrix Fails Bounds Check", __FILE__, __FUNCTION__, __LINE__ );
  Assert( mN <= (m.mN+j), "Inserted Matrix Fails Bounds Check", __FILE__, __FUNCTION__, __LINE__ );
  Assert( i >= 0, "Negative Insertion Index", __FILE__, __FUNCTION__, __LINE__ );
  Assert( j >= 0, "Negative Insertion Index", __FILE__, __FUNCTION__, __LINE__ );
  
  for ( int I = 0; I < m.mM; ++I )
    for ( int J = 0; J < m.mN; ++J )
      mData[Index(i+I,j+J)] = m(I,J);
}

void LinearAlgebra::Matrix::StraightCopy( const SparseMatrix& m )
{
  SetValues(0.0);
  for ( int i = 0; i < m.mN; ++i )
  {
    SparseMatrix::RowIterator End = m.mRows[i].end();
    for ( SparseMatrix::RowIterator iter = m.mRows[i].begin(); iter != End; ++iter )
      mData[ Index(i, iter->first) ] = iter->second;
  } 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief SparseMatrix-Matrix multiplication operator.
/// \warning The dense matrix must be square and the same size as the sparse matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Matrix LinearAlgebra::operator* ( const LinearAlgebra::SparseMatrix& A, const LinearAlgebra::Matrix& B )
{
  Assert( B.SizesMatch(A.mN, A.mN), "Sizes do not match", __FILE__, __FUNCTION__, __LINE__  );
  LinearAlgebra::Matrix M( B.mM, B.mN );
  M.SetValues(0.0);
  
  for ( int i = 0; i < M.mN; ++i )
    for ( int j = 0; j < M.mN; ++j )
    {
      SparseMatrix::RowIterator End = A.mRows[i].end();
      for ( SparseMatrix::RowIterator iter = A.mRows[i].begin(); iter != End; ++iter )
        M(i,j) += (iter->second)*B(iter->first,j);
    }
      
  return( M );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Matrix-Vector multiplication operator.
/// \warning The number of columns of the matrix must match then number length of the vector.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Vector LinearAlgebra::operator* ( const LinearAlgebra::Matrix& A, const LinearAlgebra::Vector& x )
{
  Assert( (A.mN == x.mSize), "Sizes do not match.", __FILE__, __FUNCTION__, __LINE__  );
  LinearAlgebra::Vector result( x.mSize );
  
  for ( int i = 0; i < A.mM; ++i )
    result.mData[i] = std::inner_product( A.mData + A.Index(i,0), A.mData + A.Index(i+1,0), x.mData, 0.0 );
  
  return( result );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the trace of the matrix.
///
/// Computes the trace 
/// \f[ \text{tr} \left( M \right ) = \sum_i M_{i,i} \f]
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
double LinearAlgebra::trace( const LinearAlgebra::Matrix& M )
{
  Assert( M.IsSquare(), "Tried to Calculate the Trace of a Non-Square Matrix", __FILE__, __FUNCTION__, __LINE__  );
  double Sum = 0.0;
  for ( int i = 0; i < M.mN; ++i )
    Sum += M(i,i);
  return( Sum );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the transpose of the matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Matrix LinearAlgebra::transpose( const Matrix& M )
{
  LinearAlgebra::Matrix result( M.mN, M.mM );
  for ( int i = 0; i < M.mN; ++i )
    for ( int j = 0; j < M.mM; ++j )
      result(i,j) = M(j,i);
  return( result );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Changes the current values of the matrix to the values of a matrix with a known inverse and eigenvalues.
///
/// The n-th row and the n-th column of the inverse are the set 1, 2, 3, . . . , n.  The
/// eigenvalues are all one except two, which are \f$ 6/( p (n+1) )\f$ and \f$ p/( n (5-2n) ) \f$
/// where \f$ n \f$ is the size of the matrix and \f$ p = 3 + \sqrt ( 3 (4n - 3) (n-1)/(n+1) ) \f$.
/// \warning This will not work unless the matrix is square.
/// \todo Test this function.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::Matrix::MakeTestMatrix( void )
{
  double c = mN*(mN+1.0)*(2.0*mN-5.0)/6.0;
  double d = 1.0/c;
  
  mData[ Index(mN-1, mN-1) ] = -d;
  
  for ( int i = 0; i < (mN-1); ++i )
  {
    mData[ Index(i, mN-1) ] = mData[ Index(mN-1, i) ] = (i+1)*d;
    mData[ Index(i, i) ] = d * (c - (i+1)*(i+1) );
    for ( int j = 0; j < (i - 1); ++j )
      mData[ Index( i, j ) ] = -d * (i+1) * (j+1);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the outer product of two vectors.
///
/// Computes an MxN matrix \f$ A_{i,j} = u_i * v_j \f$
/// \param u [in] Vector of size M.
/// \param v [in] Vector of size N.
/// \return An MxN matrix that is the outer product of u and v.
/// \todo Test this function.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Matrix LinearAlgebra::outerProduct( const LinearAlgebra::Vector& u, const LinearAlgebra::Vector& v )
{
  const size_t M = u.GetSize();
  const size_t N = v.GetSize();
  
  LinearAlgebra::Matrix A( M, N );
  
  for ( int i = 0; i < M; ++i )
    for ( int j = 0; j < N; ++j )
      A(i,j) = u(i) * v(j);
  
  return( A );
}

LinearAlgebra::Matrix LinearAlgebra::exp( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdExp(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::sqrt( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdSqrt(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::cbrt( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdCbrt(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::pow( const Matrix& M, const double c )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdPowx(int(M.mSize), M.mData, c, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::ln( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdLn(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::log10( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdLog10(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::cos( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdCos(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::sin( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdSin(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

void LinearAlgebra::sincos( const Matrix& M, Matrix& C, Matrix& S )
{
  #ifdef _USE_INTEL_VML_
    Assert(M.mSize == C.mSize, "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
    Assert(M.mSize == S.mSize, "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );      
    vdSinCos(int(M.mSize), M.mData, S.mData, C.mData);  
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::tan( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdTan(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::acos( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdAcos(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::asin( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdAsin(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::atan( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdAtan(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::cosh( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdCosh(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::sinh( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdSinh(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::tanh( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdTanh(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::acosh( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdAcosh(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::asinh( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdAsinh(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::atanh( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdAtanh(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::erf( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdErf(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}

LinearAlgebra::Matrix LinearAlgebra::erfc( const Matrix& M )
{
  #ifdef _USE_INTEL_VML_
    LinearAlgebra::Matrix result(M.mM, M.mN);
    vdErfc(int(M.mSize), M.mData, result.mData);
    return(result);
  #endif
}
