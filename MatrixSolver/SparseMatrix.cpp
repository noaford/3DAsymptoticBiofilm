#include "SparseMatrix.hpp"

#include "DiagonalMatrix.hpp"
#include "Lapack.hpp"
#include "LU.hpp"
#include "SafeDelete.hpp"
#include "LinearSolver.hpp"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sparse matrix inplace addition operator
/// \warning The dimension of the source and destination matrices must match
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::SparseMatrix& LinearAlgebra::SparseMatrix::operator+= ( const LinearAlgebra::SparseMatrix& m )
{
  Assert( mN == m.mN, "Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  for ( int i = 0; i < m.mN; ++i )
  {
    RowIterator End = m.mRows[i].end();
    for ( RowIterator iter = m.mRows[i].begin(); iter != End; ++iter )
      mRows[i][iter->first] += iter->second;
  }
  return( *this );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sparse matrix inplace subtraction operator
/// \warning The dimension of the source and destination matrices must match
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::SparseMatrix& LinearAlgebra::SparseMatrix::operator-= ( const LinearAlgebra::SparseMatrix& m )
{
  Assert( mN == m.mN, "Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  for ( int i = 0; i < m.mN; ++i )
  {
    RowIterator End = m.mRows[i].end();
    for ( RowIterator iter = m.mRows[i].begin(); iter != End; ++iter )
      mRows[i][iter->first] -= iter->second;
  }
  return( *this );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sparse matrix inplace right multiplication operator /f$ A = A*M /f$
/// \warning The dimension of the source and destination matrices must match
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::SparseMatrix&	LinearAlgebra::SparseMatrix::operator*= ( const LinearAlgebra::SparseMatrix& M )
{ 
  Assert( (mN == M.mN), "Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  MatrixRow* pRows = new MatrixRow[mN];
  Assert( (pRows != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  
  for ( int i = 0; i < mN; ++i )
  {
    RowIterator End = mRows[i].end();
    for ( RowIterator iter = mRows[i].begin(); iter != End; ++iter )
    {
      RowIterator bEnd = M.mRows[iter->first].end();
      for ( RowIterator biter = M.mRows[iter->first].begin(); biter != bEnd; ++biter )
        pRows[i][biter->first] += iter->second*biter->second;
    }
  }
    
  std::swap( mRows, pRows );
  SafeDeleteArray( pRows );
  
  return( *this );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Scalar inplace multiplication operator
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::SparseMatrix&	LinearAlgebra::SparseMatrix::operator*= ( const double c )
{
  if ( c == 0 )
    ClearMatrix();
  else
    for ( int i = 0; i < mN; ++i )
    {
      RowIterator End = mRows[i].end();
      for ( RowIterator iter = mRows[i].begin(); iter != End; ++iter )
        iter->second *= c;
    }
  return( *this );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Extends the number of rows/columns of the matrix.
/// 
/// The extension preserves the values within the matrix.
/// \warning Does not shrink the matrix.  If new size is less than the old size, then nothing is done.
/// \param nSize [in] The new size of the matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::Extend( const int nSize )
{
  if ( nSize > mN )
  {
    MatrixRow* temp = new MatrixRow[nSize];
    Assert( (temp != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
    
    for ( int i = 0; i < mN; ++i )
      temp[i] = mRows[i];
    
    std::swap( temp, mRows );
    SafeDeleteArray( temp );
    
    mN = nSize;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Resizes the sparse matrix
///
/// The resizing does not preserve the values within the matrix.
/// \param nSize [in] The new size of the matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::Resize( const int nSize )
{
  ReleaseMemory();
  mN = nSize;
  AllocateMemory();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Inserts a sparse matrix.
/// \todo Add range checking.
/// \param m [in] Matrix to insert.
/// \param i [in] Row to start inserting matrix into.
/// \param j [in] Column to start inserting matrix into.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::InsertBlock( const LinearAlgebra::SparseMatrix& m, const int i, const int j )
{
  for ( int I = 0; I < m.mN; ++I )
  {
    RowIterator End = m.mRows[I].end();
    for ( RowIterator iter = m.mRows[I].begin(); iter != End; ++iter )
    {
      const int C = j + iter->first;
      mRows[i+I][C] = iter->second;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Inserts a dense matrix.
/// \todo Add range checking.
/// \todo Add checking for zero elements.
/// \param m [in] Matrix to insert.
/// \param i [in] Row to start inserting matrix into.
/// \param j [in] Column to start inserting matrix into.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::InsertBlock( const LinearAlgebra::Matrix& m, const int i, const int j )
{
  for ( int I = 0; I < m.mM; ++I )
    for ( int J = 0; J < m.mN; ++J )
    {
      const int C = j + J;
      mRows[i+I][C] = m(I,J);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Returns the number of non-zero elements in the sparse matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int LinearAlgebra::SparseMatrix::GetNumNonZeroElements( void ) const
{
  int Sum = 0;
  
  for (int i = 0; i < mN; ++i )
    Sum += mRows[i].size();
  
  return( Sum );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Creates a compressed format sparse matrix for use with other libraries.
///
/// \warning The memory is owned by the caller and must be explicily released.
/// \warning The function does not release any prealloced memory.
/// \param A [out] Refernce to an array of the non-zero matrix values.  Size is equal to the # of non-zero values obtained from 
///                SparseMatrix::GetNumNonZeroElements()
/// \param ia [out] Reference to an array of the start position for each row in the A & ja arrays.  Size is equal to the # of 
///                 rows in the matrix.
/// \param ja [out] Reference to an array of the column indices for each non-zero value in the A array.  Size is equal to the # 
///                 of non-zero values obtained from SparseMatrix::GetNumNonZeroElements()
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::GetCompressedFormatArrays( double*& A, int*& ia, int*& ja ) const
{
  const int nElements = GetNumNonZeroElements();
  A = new double[nElements];
  ia = new int[mN+1];
  ja = new int[nElements];
  
  Assert( (A != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  Assert( (ia != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  Assert( (ja != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  
  int I = 1;
  for ( int i = 0; i < mN; ++i )
  {
    ia[i] = I;
    
    RowIterator End = mRows[i].end();
    for ( RowIterator iter = mRows[i].begin(); iter != End; ++iter )
    {
      A[I-1] = iter->second;
      ja[I-1] = iter->first + 1;
      ++I;
    }        
  }
  
  ia[mN] = I;  
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Creates a compressed row format sparse matrix for use with other libraries.
///
/// \warning The memory is owned by the caller and must be explicily released.
/// \warning The function does not release any prealloced memory.
/// \param A [out] Refernce to an array of the non-zero matrix values.  Size is equal to the # of non-zero values obtained from 
///                SparseMatrix::GetNumNonZeroElements()
/// \param ia [out] Reference to an array of the start position for each row in the A & ja arrays.  Size is equal to the # 
///                 of non-zero values obtained from SparseMatrix::GetNumNonZeroElements()
/// \param ja [out] Reference to an array of the column indices for each non-zero value in the A array.  Size is equal to the # of 
///                 rows in the matrix.
/// \param idxOffset [in] Starting index for array (0 for C/C++, 1 for Fortran).  Default is 0.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::GetCompressedRowArrays( double*& A, int*& ia, int*& ja, const int idxOffset ) const
{
  const int nElements = GetNumNonZeroElements();
  A = new double[nElements];
  ia = new int[nElements];
  ja = new int[mN+1];
  
  Assert( (A != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  Assert( (ia != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  Assert( (ja != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  
  int J = idxOffset;
  for ( int i = 0; i < mN; ++i )
  {
    ja[i] = J;
    
    RowIterator End = mRows[i].end();
    for ( RowIterator iter = mRows[i].begin(); iter != End; ++iter )
    {
      A[J-idxOffset] = iter->second;
      ia[J-idxOffset] = iter->first + idxOffset;
      ++J;
    }        
  }
  
  ja[mN] = J;    
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Obtains a column major array of the matrix values
/// \warning The memory is owned by the caller and must be explicitly released.
/// \warning The function does not release any preallocated memory
/// \param A [out] Column major array of matrix values.  The length of the the array is \f$ N^2 \f$.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::GetDenseArray( double*& A ) const
{
  A = new double[mN*mN];
  Assert( (A != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  
  memset( A, 0, sizeof(double)*mN*mN );
  
  for ( int i = 0; i < mN; ++i )
  {
    SparseMatrix::RowIterator End = mRows[i].end();
    for ( SparseMatrix::RowIterator iter = mRows[i].begin(); iter != End; ++iter )
      A[ i + (iter->first)*mN ] = iter->second;
  }
}

void LinearAlgebra::SparseMatrix::CopyMatrix( const LinearAlgebra::Matrix& m )
{
  Assert( (m.mM == m.mN), "Trying to Copy a Non-Square Matrix", __FILE__, __FUNCTION__, __LINE__ );
  for ( int i = 0; i < m.mM; ++i )
    for ( int j = 0; j < m.mN; ++j )
      if ( !IsZero( m(i,j) ) )
        mRows[i][j] = m(i,j);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Sparse matrix output stream insert operator
///
/// Inserts a sparse matrix into an output stream with the following format:
/// 
/// Row Index (tab) Column Index (tab) Element Value
/// \param o [in, out] Reference to the output stream.
/// \param m [in] Source matrix.
/// return Reference to the output stream after insertion.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::ostream& LinearAlgebra::operator<< ( std::ostream& o, const LinearAlgebra::SparseMatrix& m )
{
  for ( int i = 0; i < m.mN; ++i )
  {
    LinearAlgebra::SparseMatrix::RowIterator End = m.mRows[i].end();
    for ( LinearAlgebra::SparseMatrix::RowIterator iter = m.mRows[i].begin(); iter != End; ++iter )
      o << i+1 << "\t" << iter->first+1 << "\t" << iter->second << std::endl;
  }
  return(o);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief SparseMatrix-Vector multiplication operator
/// \warning The length of the vector must match the size of the sparse matrix
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::Vector LinearAlgebra::operator* ( const LinearAlgebra::SparseMatrix& A, const LinearAlgebra::Vector& x )
{
  Assert( A.mN == x.mSize, "Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  LinearAlgebra::Vector b(x.mSize);  
  
  for ( int i = 0; i < x.mSize; ++i )
  {
    LinearAlgebra::SparseMatrix::RowIterator End = A.mRows[i].end();
    for ( LinearAlgebra::SparseMatrix::RowIterator iter = A.mRows[i].begin(); iter != End; ++iter )
      b(i) += (iter->second)*x(iter->first);
  }
  return( b );
}

LinearAlgebra::Vector LinearAlgebra::trans_mult( const LinearAlgebra::SparseMatrix& A, const LinearAlgebra::Vector& x )
{
  const size_t nSize = x.GetSize();
  Assert( A.mN == nSize, "Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  LinearAlgebra::Vector b(nSize);  
  
  for ( int i = 0; i < nSize; ++i )
  {
    LinearAlgebra::SparseMatrix::RowIterator End = A.mRows[i].end();
    for ( LinearAlgebra::SparseMatrix::RowIterator iter = A.mRows[i].begin(); iter != End; ++iter )
      b(iter->first) += (iter->second)*x(i);
  }
  return( b );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the Frobenius norm of the matrix
///
/// Computes
/// \f[ \left| \left| M \right| \right| = \sqrt({\sum_i \sum_j M_{i,j}^2 } \f]
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double LinearAlgebra::frobnorm( const LinearAlgebra::SparseMatrix& M )
{
  double norm = 0.0;
  for ( int i = 0; i < M.mN; ++i )
  {
    LinearAlgebra::SparseMatrix::RowIterator End = M.mRows[i].end();
    for ( LinearAlgebra::SparseMatrix::RowIterator iter = M.mRows[i].begin(); iter != End; ++iter )
      norm += (iter->second)*(iter->second);
  }
  return( ::sqrt(norm) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the \f$ L_{\inf} \f$ norm of the matrix
///
/// Computes
/// \f[ \left| \left| M \right| \right|_{\inf} = max_i \sum_j M_{i,j} \f]
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double LinearAlgebra::infnorm( const LinearAlgebra::SparseMatrix& M )
{
  double MaxCol = 0.0;
  
  for ( int i = 0; i < M.mN; ++i )
  {
    double ColSum = 0.0;
    for ( int j = 0; j < M.mN; ++j )
      ColSum += fabs( M(i,j) );
    MaxCol = std::max(ColSum, MaxCol); 
  }
  return( MaxCol );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the \f$ L_1 \f$ norm of the matrix
/// 
/// Computes 
/// \f[ \left| \left| M \right| \right|_1 = max_j \sum_i M_{i,j} \f]
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double LinearAlgebra::onenorm( const LinearAlgebra::SparseMatrix& M )
{
  double MaxRow = 0.0; 
  for ( int i = 0; i < M.mN; ++i )
  {
    double RowSum = 0.0;
    LinearAlgebra::SparseMatrix::RowIterator End = M.mRows[i].end();
    for ( LinearAlgebra::SparseMatrix::RowIterator iter = M.mRows[i].begin(); iter != End; ++iter )
      RowSum += fabs(iter->second);
    MaxRow = std::max(RowSum, MaxRow);
  }
  return( MaxRow );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the transpose of the matrix.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::SparseMatrix LinearAlgebra::transpose( const LinearAlgebra::SparseMatrix& M )
{
  LinearAlgebra::SparseMatrix T(M.mN);
  for ( int i = 0; i < M.mN; ++i )
  {
    LinearAlgebra::SparseMatrix::RowIterator End = M.mRows[i].end();
    for ( LinearAlgebra::SparseMatrix::RowIterator iter = M.mRows[i].begin(); iter != End; ++iter )
      T.mRows[iter->first][i] = iter->second;
  }
  return( T );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Estimates the condition number of the sparse matrix.
/// \param M [in] Source matrix.
/// \param Norm [in] Character that determines which norm to use.  'I' for the \f$L_{\inf}\f$ norm, else it uses the \f$L_1\f$ norm.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double LinearAlgebra::cond( const LinearAlgebra::SparseMatrix& M, char Norm )
{
  double* A = NULL;
  M.GetDenseArray(A);               // Get Dense Matrix Version in FORTAN format for LAPACK
  double norm = 0.0;                // The norm of A
  double rcond = 0.0;               // Reciprical of the Condition Number
  int pivot[M.mN];                  // Pivot array for A = P*L*U
  double DoubleWorkArray[4*M.mN];   // Workspace array of doubles for LAPACK
  int IntWorkArray[M.mN];           // Workspace array of ints for LAPACK
  int Info = 0;                     // Reports if the function worked.
  int nSize = int(M.mN);                 // Matrix size
  
  if ( Norm == 'I' )
    norm = infnorm(M);              // Compute the infinity norm of A
  else
  {
    norm = '1';
    norm = onenorm(M);              // Compute the 1-norm of A
  }
  
  dgetrf( &nSize, &nSize, A, &nSize, pivot, &Info );
  
  if ( Info < 0)
  {
    switch( Info )
    {
    case (-1):
      std::cerr << "\t***** Error in dgecon - Bad norm type." << std::endl;
      break;
      
    case (-2):
    std::cerr << "\t***** Error in dgecon - Invalid matrix size." << std::endl;
      break;
      
    case (-3):
    std::cerr << "\t***** Error in dgecon - Invalid matrix array." << std::endl;
      break;
      
    case (-4):
    std::cerr << "\t***** Error in dgecon - Invalid dimension size." << std::endl;
      break;
      
    case (-5):
    std::cerr << "\t***** Error in dgecon - Invalid matrix norm." << std::endl;
      break;
      
    case (-6):
    std::cerr << "\t***** Error in dgecon - Bad rcond parameter." << std::endl;
      break;
      
    case (-7):
    std::cerr << "\t***** Error in dgecon - Invalid double workspace." << std::endl;
      break;
      
    case (-8):
    std::cerr << "\t***** Error in dgecon - Invalid integer workspace." << std::endl;
      break;
      
    case (-9):
      std::cerr << "\t***** Error in dgecon - Bad info parameter." << std::endl;
      break;
      
    default:
      std::cerr << "\t***** Unknown error in dgecon." << std::endl;  
    };
  }
  else
    if ( Info > 0 )
      std::cerr << "\t***** Warnging: In dgetrf, U(" << Info << ", " << Info << ") is exactly zero." << std::endl;
  
  dgecon( &Norm, &nSize, A, &nSize, &norm, &rcond, DoubleWorkArray, IntWorkArray, &Info );
  
  if ( Info )
  {
    switch( Info )
    {
    case (-1):
      std::cerr << "\t***** Error in dgecon - Bad norm type." << std::endl;
      break;
      
    case (-2):
    std::cerr << "\t***** Error in dgecon - Invalid matrix size." << std::endl;
      break;
      
    case (-3):
    std::cerr << "\t***** Error in dgecon - Invalid matrix array." << std::endl;
      break;
      
    case (-4):
    std::cerr << "\t***** Error in dgecon - Invalid dimension size." << std::endl;
      break;
      
    case (-5):
    std::cerr << "\t***** Error in dgecon - Invalid matrix norm." << std::endl;
      break;
      
    case (-6):
    std::cerr << "\t***** Error in dgecon - Bad rcond parameter." << std::endl;
      break;
      
    case (-7):
    std::cerr << "\t***** Error in dgecon - Invalid double workspace." << std::endl;
      break;
      
    case (-8):
    std::cerr << "\t***** Error in dgecon - Invalid integer workspace." << std::endl;
      break;
      
    case (-9):
      std::cerr << "\t***** Error in dgecon - Bad info parameter." << std::endl;
      break;
      
    default:
      std::cerr << "\t***** Unknown error in dgecon." << std::endl;  
    };
  }
  
  SafeDeleteArray(A);
  
  return( rcond );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Scales a column in the sparse matrix by a given scaling factor.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::ScaleColumn( const int j, const double s )
{
  for ( int i = 0; i < mN; ++i )
  { 
    RowIterator iter = mRows[i].find(j); 
    if ( iter != mRows[i].end() ) 
      iter->second *= s; 
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes scaling factors to reduce the condition number of the sparse matrix
/// \param rows [out] Vector of row scaling factors.
/// \param cols [out] Vector of column scaling factors.
/// \param rowRatio [out] The ratio \f$ \frac{ \min_i \max_j M_{i,j} }{ \max_i \max_j M_{i,j} } \f$
/// \param colRatio [out] The ratio \f$ \frac{ \min_j \max_i M_{i,j} }{ \max_j \max_i M_{i,j} } \f$
/// \param Max [out] Absolute maximum element value.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::ComputeScalingFactors( Vector& rows, Vector& cols, double& rowRatio, double& colRatio, double& Max ) const
{
  rows = Vector(mN, 0.0);
  cols = Vector(mN, 0.0);
  Max = 0.0;  
  
  if ( mN == 0 )
    return;
  
  for ( int i = 0; i < mN; ++i )
  {
    RowIterator RowEnd = mRows[i].end();
    for ( RowIterator iter = mRows[i].begin(); iter != RowEnd; ++iter )
    {
      const double Value = fabs( iter->second );
      rows(i) = std::max( rows(i), Value );
    }
  }
  
  const double RowMax = absmax(rows);
  const double RowMin = absmin(rows);
  
  Max = RowMax;
  
  Assert( RowMax != 0.0, "Zero Row in Matrix", __FILE__, __FUNCTION__, __LINE__ );

  for ( int i = 0; i < mN; ++i )
    rows(i) = 1.0 / rows(i);
  
  rowRatio = RowMin/RowMax;
  
  for ( int j = 0; j < mN; ++j )
    for ( int i = 0; i < mN; ++i )
    {
      const double Value = fabs( mRows[i][j] ) * rows(i);
      cols(j) = std::max( cols(j), Value  );
    }

  const double ColMax = max(cols);
  const double ColMin = min(cols);
  
  Assert( ColMax != 0.0, "Zero Column in Matrix", __FILE__, __FUNCTION__, __LINE__ );
  
  for ( int j = 0; j < mN; ++j )
    cols(j) = 1.0/cols(j);
  
  colRatio = ColMin/ColMax;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Takes the dot product of a row in a sparse matrix with a vector
/// \warning the length of the vector must match the number of columns in the matrix.
/// \param M [in] Source matrix.
/// \param v [in] Source vector.
/// \param Index [in] Row index in the sparse matrix.
/// \return \f$ \sum_j M_{i,j} * v_j \f$
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double LinearAlgebra::RowMult( const LinearAlgebra::SparseMatrix& M, const LinearAlgebra::Vector& v, const int Index )
{
  double x = 0.0;
  LinearAlgebra::SparseMatrix::RowIterator End = M.mRows[Index].end();
  for ( LinearAlgebra::SparseMatrix::RowIterator iter = M.mRows[Index].begin(); iter != End; ++iter )
    x += (iter->second)*v(iter->first);

  return( x );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Creates a coordinate format sparse matrix for use with other libraries.
///
/// \warning The memory is owned by the caller and must be explicily released.
/// \warning The function does not release any prealloced memory.
/// \param A [out] Refernce to an array of the non-zero matrix values.  Size is equal to the # of non-zero values obtained from 
///                SparseMatrix::GetNumNonZeroElements()
/// \param ia [out] Reference to an array of the row index for each element in A. Size is equal to the # of non-zero values 
///                 obtained from SparseMatrix::GetNumNonZeroElements()
/// \param ja [out] Reference to an array of the column index for each element in A. Size is equal to the # of non-zero values 
///                 obtained from SparseMatrix::GetNumNonZeroElements()
/// \param idxOffset [in] Starting index for array (0 for C/C++, 1 for Fortran).  Default is 0.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LinearAlgebra::SparseMatrix::GetCoordinateFormatArrays( double*& A, int*& ia, int*& ja, const int idxOffset ) const
{
  const int nElements = GetNumNonZeroElements();
  A  = new double[nElements];
  ia = new int[nElements];
  ja = new int[nElements];
  
  Assert( (A != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  Assert( (ia != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  Assert( (ja != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
    
  int I = 0;
  for ( int i = 0; i < mN; ++i )
  {        
    RowIterator End = mRows[i].end();
    for ( RowIterator iter = mRows[i].begin(); iter != End; ++iter )
    {
       A[I] = iter->second;
      ja[I] = iter->first + idxOffset;
      ia[I] = i + idxOffset; 
      ++I;
    }        
  }    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Estimates the dominate eigenvalue of the sparse matrix using the power method.
/// \param M [in] Source matrix.
/// \param x [in, out] Initial guess. Hold an estimate of the dominate eigenvector after the function is called.  
///                    Should not be zero or iteration will not succeed.
/// \param epsilon [in] Relative convergence tolerance.
/// \param maxIter [in] Maximum iterations of the power method.
/// \return Approximation of the dominate eigenvalue of the sparse matrix.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double LinearAlgebra::maxEV( const SparseMatrix& M, Vector& x, const double epsilon, const int maxIter )
{
  Vector cX = x;                        // Current x
  Vector nX = M*x;                      // Next x
  double cL = 0.0;                      // Current lambda
  double nL = dot(cX, nX)/mag2(cX);     // Next lambda
  
  int numIter = 1;
  double delta = fabs(nL - cL);
  double deltaX = mag(nX - cX);
  
  while ( ( delta > epsilon ) && ( deltaX > epsilon ) && ( numIter <= maxIter ) )
  {
    cX = nX;                            // Set the next x as current.
    cX /= absmax(cX);                   // Scale the current vector to keep the value reasonable.
    
    nX = M*cX;                          // Get the next x.
    
    cL = nL;                            // Set the next lambda as current.
    nL = dot(cX, nX)/mag2(cX);          // Get the next lambda.
    
    ++numIter;                          // Increment the number of iterations.
    delta = nL - cL;                    // Compute the change.
    deltaX = mag( nX - cX );            
        
  }
    
  if ( numIter > maxIter )
  {
    std::cerr << "*** Maximum iterations reached when estimating the dominate eigenvalue!\n";
    std::cerr << "***       Eigenvalue Delta = " << delta << "\n";
    std::cerr << "***      Eigenvector Delta = " << deltaX << std::endl;
  }
  
  x = nX/nL;
  
  return( nL );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Estimates the absolute smallest eigenvalue of the sparse matrix using the inverse power method.
/// \param M [in] Source matrix.
/// \param x [in, out] Initial guess. Hold an estimate of the eigenvector after the function is called.  
/// \param epsilon [in] Relative convergence tolerance.
/// \param maxIter [in] Maximum iterations of the power method.
/// \return Approximation of the absolute smalled eigenvalue of the sparse matrix.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double LinearAlgebra::minEV( const SparseMatrix& M, Vector& x, const double epsilon, const int maxIter )
{
  Vector cX = x;                        // Current x
  Vector nX(x.GetSize());               // Next x
  
  LU luM( M );
  
  //if ( !DirectSolver::Solve( M, cX, nX, false ) )  // Solve for next x
  if ( !luM.Solve( cX, nX ) )
  { // Solver failed
    std::cerr << "*** Unable to linear system to calculate the eigenvalue!\n" << std::endl;
    return( 0.0 );
  }
  
  double cL = 0.0;                      // Current lambda
  double nL = dot(cX, nX)/mag2(cX);     // Next lambda
  
  int numIter = 1;
  double delta = fabs(nL - cL);
  double deltaX = mag(nX - cX);
  
  while ( ( delta > epsilon ) && ( deltaX > epsilon ) && ( numIter <= maxIter ) )
  {
    cX = nX;                            // Set the next x as current.
    cX /= absmax(cX);                   // Scale the current vector to keep the value reasonable.
    
    //if ( !DirectSolver::Solve( M, cX, nX, false ) )  // Solve for next x
    if ( !luM.Solve( cX, nX ) )
    { // Solver failed
      std::cerr << "*** Unable to linear system to calculate the eigenvalue!\n" << std::endl;
      return( 0.0 );
    }
    
    cL = nL;                            // Set the next lambda as current.
    nL = dot(cX, nX)/mag2(cX);          // Get the next lambda.
    
    ++numIter;                          // Increment the number of iterations.
    delta = fabs( nL - cL );            // Compute the change
    deltaX = mag( nX - cX );            
  }
    
  if ( numIter > maxIter )
  {
    std::cerr << "*** Maximum iterations reached when estimating the absolute smallest eigenvalue!\n";
    std::cerr << "***       Eigenvalue Delta = " << delta << "\n";
    std::cerr << "***      Eigenvector Delta = " << deltaX << std::endl;
  }
  
  x = nX/nL;
  
  return( 1.0/nL );
}

double LinearAlgebra::SparseMatrix::GetMaxElement( void ) const
{
  double max = GetMaxRowElement(0);
  
  for ( int i = 1; i < mN; ++i )
    max = std::max( max, GetMaxRowElement(i) );
    
  return( max );  
}

double LinearAlgebra::SparseMatrix::GetMaxRowElement( const int i, const int offset ) const
{
  Assert( i < mN, "Invalid Row Index", __FILE__, __FUNCTION__, __LINE__ );
  Assert( i >= 0, "Invalid Row Index", __FILE__, __FUNCTION__, __LINE__ ); 
  Assert( offset < mN, "Invalid Offset", __FILE__, __FUNCTION__, __LINE__ );
  Assert( offset >= 0, "Invalid Offset", __FILE__, __FUNCTION__, __LINE__ );
  
  return( std::max_element( mRows[i].begin(), mRows[i].end() )->second );
}

double LinearAlgebra::SparseMatrix::GetMaxColumnElement( const int j, const int offset ) const
{
  Assert( j < mN, "Invalid Column Index", __FILE__, __FUNCTION__, __LINE__ ); 
  Assert( j >= 0, "Invalid Column Index", __FILE__, __FUNCTION__, __LINE__ ); 
  Assert( offset < mN, "Invalid Offset", __FILE__, __FUNCTION__, __LINE__ );
  Assert( offset >= 0, "Invalid Offset", __FILE__, __FUNCTION__, __LINE__ );
  
  double max = -std::numeric_limits<double>::max();
  
  for ( int i = 0; i < mN; ++i )
    max = std::max( max, (*this)(i,j) );
    
  return( max );
}

double LinearAlgebra::SparseMatrix::GetAbsMaxElement( void ) const
{
  double max = GetAbsMaxRowElement(0);
  
  for ( int i = 1; i < mN; ++i )
    max = std::max( max, GetAbsMaxRowElement(i) );
    
  return( max );  
}

double LinearAlgebra::SparseMatrix::GetAbsMaxRowElement( const int i, const int offset ) const
{
  Assert( i < mN, "Invalid Row Index", __FILE__, __FUNCTION__, __LINE__ );
  Assert( i >= 0, "Invalid Row Index", __FILE__, __FUNCTION__, __LINE__ ); 
  Assert( offset < mN, "Invalid Offset", __FILE__, __FUNCTION__, __LINE__ );
  Assert( offset >= 0, "Invalid Offset", __FILE__, __FUNCTION__, __LINE__ );
  
  double max = 0.0;
  
  const RowIterator RowEnd = mRows[i].end(); 
  
  for ( RowIterator iter = mRows[i].begin(); iter != RowEnd; ++iter )
    max = std::max( max, fabs(iter->second) );
  
  return( max );
}

double LinearAlgebra::SparseMatrix::GetAbsMaxColumnElement( const int j, const int offset ) const
{
  Assert( j < mN, "Invalid Column Index", __FILE__, __FUNCTION__, __LINE__ ); 
  Assert( j >= 0, "Invalid Column Index", __FILE__, __FUNCTION__, __LINE__ ); 
  Assert( offset < mN, "Invalid Offset", __FILE__, __FUNCTION__, __LINE__ );
  Assert( offset >= 0, "Invalid Offset", __FILE__, __FUNCTION__, __LINE__ );
  
  double max = 0.0;
  
  for ( int i = offset; i < mN; ++i )
    max = std::max( max, fabs( (*this)(i,j) ) );
    
  return( max );
}

LinearAlgebra::SparseMatrix& LinearAlgebra::SparseMatrix::operator+= ( const LinearAlgebra::DiagonalMatrix& m )
{
  Assert( SizesMatch( m.GetSize() ), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  for ( int i = 0; i < mN; ++i )
    mRows[i][i] += m(i);
    
  return( *this );
}

LinearAlgebra::SparseMatrix& LinearAlgebra::SparseMatrix::operator-= ( const LinearAlgebra::DiagonalMatrix& m )
{
  Assert( SizesMatch( m.GetSize() ), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  for ( int i = 0; i < mN; ++i )
    mRows[i][i] -= m(i);
    
  return( *this );
}

LinearAlgebra::SparseMatrix& LinearAlgebra::SparseMatrix::operator*= ( const LinearAlgebra::DiagonalMatrix& m )
{
  Assert( SizesMatch( m.GetSize() ), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  for ( int i = 0; i < mN; ++i )
    mRows[i][i] *= m(i);
    
  return( *this );
}

LinearAlgebra::SparseMatrix LinearAlgebra::operator+ ( const LinearAlgebra::SparseMatrix& A, const LinearAlgebra::DiagonalMatrix& B )
{
  SparseMatrix result( A );
  return( result += B );
}

LinearAlgebra::SparseMatrix LinearAlgebra::operator+ ( const LinearAlgebra::DiagonalMatrix& A, const LinearAlgebra::SparseMatrix& B )
{
  SparseMatrix result( B );
  return( result += A );
}

LinearAlgebra::SparseMatrix LinearAlgebra::operator- ( const LinearAlgebra::SparseMatrix& A, const LinearAlgebra::DiagonalMatrix& B )
{
  SparseMatrix result( A );
  return( result -= B );
}

LinearAlgebra::SparseMatrix LinearAlgebra::operator- ( const LinearAlgebra::DiagonalMatrix& A, const LinearAlgebra::SparseMatrix& B )
{
  SparseMatrix result( B );
  return( result -= A );
}

LinearAlgebra::SparseMatrix LinearAlgebra::operator* ( const LinearAlgebra::SparseMatrix& A, const LinearAlgebra::DiagonalMatrix& B )
{
  SparseMatrix result( A );
  return( result *= B );
}

LinearAlgebra::SparseMatrix LinearAlgebra::operator* ( const LinearAlgebra::DiagonalMatrix& A, const LinearAlgebra::SparseMatrix& B )
{
  SparseMatrix result( B );
  return( result.LeftMultiply( A ) );
}

LinearAlgebra::SparseMatrix& LinearAlgebra::SparseMatrix::LeftMultiply( const LinearAlgebra::DiagonalMatrix& A )
{
  Assert( SizesMatch( A.GetSize() ), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );

  for ( int i = 0; i < mN; ++i )
  {
    const RowIterator RowEnd = mRows[i].end(); 
  
    for ( RowIterator iter = mRows[i].begin(); iter != RowEnd; ++iter )
      iter->second *= A(iter->first);
  }
  
  return( *this );
}

LinearAlgebra::SparseMatrix& LinearAlgebra::SparseMatrix::LeftMultiply( const LinearAlgebra::SparseMatrix& B )
{
  Assert( SizesMatch( B.mN ), "Matrix Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );  
  
  MatrixRow* pRows = new MatrixRow[mN];
  Assert( (pRows != NULL), "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  
  for ( int i = 0; i < mN; ++i )
  {
    
    RowIterator bEnd = B.mRows[i].end();
    for ( RowIterator biter = B.mRows[i].begin(); biter != bEnd; ++biter )    
    {
      RowIterator End = mRows[biter->first].end();
      for ( RowIterator iter = mRows[biter->first].begin(); iter != End; ++iter )
        pRows[i][iter->first] += iter->second*biter->second;
    }
  }
    
  std::swap( mRows, pRows );
  SafeDeleteArray( pRows );
  
  return( *this );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Computes the sum of a matrix row.
/// \param i [in] Row of the matrix to sum.
/// \return The sum of the i-th row of the matrix.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double LinearAlgebra::SparseMatrix::ComputeRowSum( const int i ) const
{
  double Sum = 0.0;
  
  RowIterator End = mRows[i].end();
  for ( RowIterator iter = mRows[i].begin(); iter != End; ++iter )
  {
    Sum += iter->second;
  }
  
  return( Sum );
}
