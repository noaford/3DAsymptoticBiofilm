#ifndef _SPARSEMATRIX_HPP_10142004_
#define _SPARSEMATRIX_HPP_10142004_

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
/// \file
/// \brief Square sparse matrix class
/// \ingroup linalg 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <map>

#include "Assert.hpp"
#include "LinearAlgebra.hpp"
#include "Matrix.hpp"

namespace LinearAlgebra
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Square sparse matrix class
  /// \ingroup linalg
  /// \todo Add range checking
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
  class SparseMatrix
  {       
    friend class LinearAlgebra::Matrix;
    friend class LinearAlgebra::SparseMatrixBlock;
    
  public:
  //@{
    /// \name Constructors & Destructor
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Default constructor.
    /// 
    /// Constructs a 1x1 zero sparse matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    SparseMatrix( void ) : mN(1), mRows(NULL) 
    { AllocateMemory(); }	
  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Size constructor
    /// 
    /// Creates a NxN zero sparse matrix.
    /// \param N [in] Number of rows and columns of the matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    explicit SparseMatrix( const size_t N ) : mN(N), mRows(NULL) 
    { AllocateMemory(); }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sparse matrix copy constructor.
    ///
    /// \param m [in] Source matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	   
    SparseMatrix( const SparseMatrix& m ) : mN(m.mN), mRows(NULL) 
    { 
      AllocateMemory(); 
      CopyMatrix(m); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Dense matrix copy constructor.
    ///
    /// Creates a sparse matrix from the non-zero values of the source matrix.
    /// \warning Source matrix must be square.
    /// \param m [in] Source matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    SparseMatrix( const Matrix& m ) : mN(m.mN), mRows(NULL)
    {
      AllocateMemory();
      CopyMatrix(m);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Destructor
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ~SparseMatrix( void ) 
    { ReleaseMemory(); }
  //@}    

  //@{
    /// \name Assignment operators.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Dense matrix assignment operator
    /// 
    /// Create a sparse matrix from the non-zero elements of the source matrix.
    /// 
    /// \warning Source matrix must be square.
    /// \param m [in] Source matrix.
    /// \return Reference to the destination matrix after the assignment.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    SparseMatrix& operator= ( const Matrix& m ) 
    { 
      mN = m.mN; 
      AllocateMemory(); 
      CopyMatrix(m); 
      return(*this); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sparse matrix assignment operator.
    /// \param m [in] Source matrix
    /// \return Reference to the destination matrix after the assignment.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////           
    SparseMatrix& operator= ( const SparseMatrix& m ) 
    { 
      mN = m.mN; 
      AllocateMemory(); 
      CopyMatrix(m); 
      return(*this); 
    }
        
    SparseMatrix& operator+= ( const SparseMatrix& m );	    
   	SparseMatrix& operator-= ( const SparseMatrix& m );            
    SparseMatrix&	operator*= ( const SparseMatrix& m );          
    SparseMatrix&	operator*= ( const double c );
    
    SparseMatrix& LeftMultiply( const SparseMatrix& m );
    
    SparseMatrix& operator+= ( const DiagonalMatrix& m );
    SparseMatrix& operator-= ( const DiagonalMatrix& m );    
    SparseMatrix& operator*= ( const DiagonalMatrix& m );
    
    SparseMatrix& LeftMultiply( const DiagonalMatrix& m );
    
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Scalar inplace division operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	          
    SparseMatrix& operator/= ( const double c )
    { return( (*this) *= (1.0/c) ); }
  //@}
  
  //@{
    /// \name Accessors
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief %Matrix size accessor
    /// \return The number of rows/columns of the matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	  
    size_t GetSize( void ) const { return( mN ); }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief %Matrix element reference operator.
    /// \param i [in] Row index.
    /// \param j [in] Column index.
    /// \return Reference to the matrix element at (i,j).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////  	 
    double& operator() ( const int i, const int j ) 
    { 
      Assert( i < mN, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      Assert( j < mN, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      Assert( i > -1, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      Assert( j > -1, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      return( mRows[i][j] ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief %Matrix element value operator.
    /// \param i [in] Row index.
    /// \param j [in] Column index.
    /// \return Value of the matrix element at (i,j).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    double operator() (const int i, const int j ) const 
    { 
      Assert( i < mN, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      Assert( j < mN, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      Assert( i > -1, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      Assert( j > -1, "Invalid Index", __FILE__, __FUNCTION__, __LINE__ );
      
      RowIterator iter = mRows[i].find(j); 
      if ( iter == mRows[i].end() ) 
        return( 0.0 ); 
      else 
        return( iter->second );
    }
    
    bool IsNonZero( const int i, const int j ) const;
    bool IsZero( const int i, const int j ) const;
    
    int GetNumNonZeroElements( void ) const;
    
    double ComputeRowSum( const int i ) const;
    
    double GetMaxElement( void ) const;
    double GetMaxRowElement( const int i, const int offset = 0 ) const;
    double GetMaxColumnElement( const int j, const int offset = 0 ) const;
    
    double GetAbsMaxElement( void ) const;
    double GetAbsMaxRowElement( const int i, const int offset = 0 ) const;
    double GetAbsMaxColumnElement( const int j, const int offset = 0 ) const;
  //@}
     
  //@{
    /// \name Functions for resizing the matrix.
    void Extend( const int nSize );       
    void Resize( const int nSize );
  //@}
    
  //@{
    /// \name Block assignment functions.
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sets a row of the matrix to a given vector.
    ///
    /// \todo Add index checking.
    /// \todo Add code to check for zeros in the vector before they are copied.
    /// \warning Copies all values including zeros.
    /// \param i [in] Destination row.
    /// \param v [in] Source vector.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    void SetRow( const int i, const Vector& v ) 
    { 
      for ( int j = 0; j < v.mSize; ++j ) 
        mRows[i][j] = v(j); 
    }
	
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sets a column of the matrix to a given vector.
    ///
    /// \todo Add index checking.
    /// \todo Add code to check for zeros in the vector before they are copied.
    /// \warning Copies all values including zeros.    
    /// \param j [in] Destination column.
    /// \param v [in] Source vector.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
    void SetColumn( const int j, const Vector& v ) 
    { 
      for ( int i = 0; i < v.mSize; ++i ) 
        mRows[i][j] = v(i); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Clears a row in the matrix.
    /// \param i [in] Row to clear (Set to zero).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      
    void ClearRow( const int i ) 
    { MatrixRow().swap(mRows[i]); }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Clears the entire matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void ClearMatrix( void )
    {
      for ( int i = 0; i < mN; ++i )
        ClearRow(i);
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Swaps two rows in the matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void SwapRows( const int a, const int b ) 
    { std::swap( mRows[a], mRows[b] ); }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Zeros out an element of the matrix.
    /// \param i [in] Row index.
    /// \param j [in] Column index.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void ZeroElement( const int i, const int j ) 
    { mRows[i].erase(j); }       
            
    void InsertBlock( const SparseMatrix& m, const int i, const int j );        
    void InsertBlock( const Matrix& m, const int i, const int j );
  //@}   
                               
  //@{
    /// \name Array functions
    void GetCompressedFormatArrays( double*& A, int*& ia, int*& ja ) const;
    void GetCompressedRowArrays( double*& A, int*& ia, int*& ja, const int idxOffset = 0 ) const;
    void GetDenseArray( double*& A ) const;
    void GetCoordinateFormatArrays( double*& A, int*& rows, int*& ja, const int idxOffset = 0 ) const;
  //@}
    
  //@{
    /// \name Scaling functions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Scales a row by a constant.
    /// \param i [in] Row index.
    /// \param s [in] Scaling constant.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void ScaleRow( const int i, const double s )
    {
      Assert( i > -1, "Negative Index", __FILE__, __FUNCTION__, __LINE__ );
      const RowIterator RowEnd = mRows[i].end(); 
      for ( RowIterator iter = mRows[i].begin(); iter != RowEnd; ++iter )
        iter->second *= s;
    }
    
    void ScaleColumn( const int j, const double s );    
    void ComputeScalingFactors( Vector& rows, Vector& cols, double& rowRatio, double& colRatio, double& max ) const;    	        
  //@}
  
  //@{
    /// \name Stream operators 
    friend std::ostream& operator<< ( std::ostream& o, const SparseMatrix& m );
  //@}
    
  
  //@{
    /// \name Elementry operators
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sparse matrix addition operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend SparseMatrix operator+ ( const SparseMatrix& A, const SparseMatrix& B ) 
    { 
      LinearAlgebra::SparseMatrix result(A); 
      return( result += B ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sparse matrix subtraction operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend SparseMatrix operator- ( const SparseMatrix& A, const SparseMatrix& B ) 
    { 
      LinearAlgebra::SparseMatrix result(A); 
      return( result -= B ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sparse matrix multiplication operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend SparseMatrix operator* ( const SparseMatrix& A, const SparseMatrix& B ) 
    {
      LinearAlgebra::SparseMatrix M(A);
      return( M *= B );
    }
    
    friend Matrix operator* ( const SparseMatrix& A, const Matrix& B );
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief SparseMatrix-Scalar mulitplication operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend SparseMatrix operator* ( const SparseMatrix& A, const double c ) 
    { return( c*A ); }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Scalar-SparseMatrix mulitplication operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend SparseMatrix operator* ( const double c, const SparseMatrix& A )
    {
      LinearAlgebra::SparseMatrix result(A);
      return( result *= c );
    }
    
    friend SparseMatrix operator+ ( const SparseMatrix& A, const DiagonalMatrix& B );
    friend SparseMatrix operator+ ( const DiagonalMatrix& A, const SparseMatrix& B );
    friend SparseMatrix operator- ( const SparseMatrix& A, const DiagonalMatrix& B );
    friend SparseMatrix operator- ( const DiagonalMatrix& A, const SparseMatrix& B );
    friend SparseMatrix operator* ( const SparseMatrix& A, const DiagonalMatrix& B );
    friend SparseMatrix operator* ( const DiagonalMatrix& A, const SparseMatrix& B );
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sparse matrix negation operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend SparseMatrix operator- ( const SparseMatrix& A ) 
    { return( -1.0*A ); }
       
    friend Vector operator* ( const SparseMatrix& A, const Vector& x );
    friend Vector trans_mult( const SparseMatrix& A, const Vector& x );
  //@}
      
  //@{
    /// \name Matrix norms
    friend double frobnorm( const SparseMatrix& M );
    
    friend double onenorm( const SparseMatrix& M );
    friend double infnorm( const SparseMatrix& M );
  //@}
    
  //@{
    /// \name Functions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Computes the trace of the matrix.
    ///
    /// Computes the trace 
    /// \f[ \text{tr} \left( M \right ) = \sum_i M_{i,i} \f]
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend double trace( const SparseMatrix& M ) 
    { 
      double Trace = 0.0; 
      for ( int i = 0; i < M.mN; ++i ) 
        Trace += M(i,i); 
      return( Trace ); 
    }
        
    friend SparseMatrix transpose( const SparseMatrix& M );
    friend double cond( const SparseMatrix& A, char Norm);
    friend double RowMult( const SparseMatrix& A, const Vector& v, const int Index );
  //@}
    
  //@{ Indentity matrix functions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Creates a NxN identity matrix.
    /// \param N [in] Size of the matrix.
    /// \return NxN identity matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend SparseMatrix identity( const int N ) 
    { 
      LinearAlgebra::SparseMatrix M( N ); 
      for ( int i = 0; i < N; ++i ) 
        M(i,i) = 1.0;  
      return( M ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Converts a matrix into an identity matrix.
    /// \param M [in, out] Source matrix
    /// \return NxN identity matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend SparseMatrix& identity( SparseMatrix& M ) 
    { 
      for ( int i = 0; i < M.mN; ++i ) 
      { 
        M.ClearRow(i); 
        M(i,i) = 1.0; 
      } 
      return( M ); 
    }
    
  //@}
  
  //@{
    // \name Eigenvalue functions    
    friend double maxEV( const SparseMatrix& M, Vector& x, const double epsilon, const int maxIter );
    friend double minEV( const SparseMatrix& M, Vector& x, const double epsilon, const int maxIter );
  //@}    

  private:    
    typedef std::pair<int, double> MatrixElement;
    typedef std::map<int, double>  MatrixRow;    
    typedef MatrixRow::iterator    RowIterator;
    
    bool SizesMatch( const size_t nSize )
    { return( nSize == mN ); }
    
    void AllocateMemory( void ) 
    { 
      if (mRows == NULL) 
        mRows = new MatrixRow[mN]; 
    }
    
    void ReleaseMemory( void ) 
    { 
      if (mRows != NULL) 
      { 
        delete[] mRows; 
        mRows = NULL; 
      } 
    }
    
    void CopyMatrix( const SparseMatrix& m ) 
    { 
      for ( int i = 0; i < m.mN; ++i )
        mRows[i] = m.mRows[i];
    }
    
    void CopyMatrix( const Matrix& m );
    
    bool IsZero( const double x ) 
    { return( fabs(x) < ( 10.0 * std::numeric_limits<double>::max() ) ); }
    
    size_t mN;
    MatrixRow* mRows;
  };   
}

#endif

