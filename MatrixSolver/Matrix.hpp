#ifndef _MATRIX_HPP_10142004_
#define _MATRIX_HPP_10142004_

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Dense matrix class.
/// \ingroup linalg
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
#include <numeric> 
#include <ostream> 
 #include <iostream>
#include "LinearAlgebra.hpp"
#include "Math.hpp"
#include "SafeDelete.hpp"
#include "Vector.hpp"

namespace LinearAlgebra
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// \brief Dense %Matrix class.
  /// \ingroup linalg
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class Matrix
  {
  friend class SparseMatrix;
    
  public:  
  /// \name Constructors & Destructors
  //@{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Default constructor.
    ///
    /// Creates a 1x1 matrix whose value is undefined.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    Matrix( void ) : mM(0), mN(0), mSize(0), mData(NULL) 
    { AllocateMemory(); }
	
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  /// \brief Matrix size constructor.
	  ///
    /// Creates a MxN matrix whose elements are undefined.
	  /// \param M [in] # of rows in the matrix.
    /// \param N [in] # of columns in the matrix.
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Matrix( const size_t M, const size_t N ) : mM(M), mN(N), mSize(M*N), mData(NULL) 
    { AllocateMemory(); }
	
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  /// \brief Initilized matrix constructor
    /// 
    /// Creates a MxN matrix whose elements are of a given value.
    /// \param M [in] # of rows in the matrix.
    /// \param N [in] # of columns in the matrix.
    /// \param Value [in] Value to assign to each matrix element.
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Matrix( const size_t M, const size_t N, double Value) : mM(M), mN(N), mSize(M*N), mData(NULL) 
    { 
      AllocateMemory(); 
      SetValues(Value); 
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Dense matrix copy constructor.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Matrix( const Matrix& m ) : mM(m.mM), mN(m.mN), mSize(m.mSize), mData(NULL) 
    { 
      AllocateMemory(); 
      StraightCopy(m); 
    }
        
    Matrix( const SparseMatrix& m ); 
     
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix destructor    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ~Matrix( void ) 
    { ReleaseMemory(); }
  //@}  
    
  /// \name Initializers
  //@{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sets all the element values in a matrix to a given value.
    /// \param Value [in] Value to assign to all the elements of the matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////        
    void SetValues( const double Value ) 
    { std::fill( mData, mData+mSize, Value ); }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Zeros out the matrix.
    ///
    /// Same as calling Matrix::SetValues(0.0).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  	void ZeroValues( void ) 
    { SetValues(0.0); }
  //@}  
    
  /// \name Assignment operators
  //@{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Dense matrix assignment operator.
    ///
    /// If the source and the destination matrix sizes match, then the value are copied, otherwise memory is reallocated in the 
    /// destination matrix to hold the values of the source matrix.
    /// \param m [in] Source matrix.
    /// \return Reference to the destination matrix after the assignment.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////           
    Matrix& operator= ( const Matrix& m ) 
    { 
      if ( SizesMatch(m.mM, m.mN) )
        StraightCopy(m);
      else 
      {
        mM = m.mM;
        mN = m.mN;
        mSize = m.mSize;
        AllocateCopy(m);
      }
      return(*this); 
    }
               
     Matrix& operator= ( const SparseMatrix& m );
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Dense matrix inplace addition operator.
    ///
    /// \warning The size of the source matrix must match the size of the destination matrix.
    /// \param m [in] Source matrix.
    /// \return Reference to the destination matrix after the addition.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Matrix&	operator+= ( const Matrix& m ) 
    { 
      Assert( m.mM == mM, "Matrix Sizes do not Match", __FILE__, __FUNCTION__, __LINE__ );
      Assert( m.mN == mN, "Matrix Sizes do not Match", __FILE__, __FUNCTION__, __LINE__ );
      std::transform( mData, mData+mSize, m.mData, mData, std::plus<double>() ); 
      return(*this);  
    }
	
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Dense matrix inplace subtraction operator.
    ///
    /// \warning The size of the source matrix must match the size of the destination matrix.
    /// \param m [in] Source matrix.
    /// \return Reference to the destination matrix after the subtraction.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   	Matrix& operator-= ( const Matrix& m ) 
    { 
      Assert( m.mM == mM, "Matrix Sizes do not Match", __FILE__, __FUNCTION__, __LINE__ );
      Assert( m.mN == mN, "Matrix Sizes do not Match", __FILE__, __FUNCTION__, __LINE__ );
      std::transform( mData, mData+mSize, m.mData, mData, std::minus<double>() ); 
      return(*this);  
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Scalar inplace multiplication operator.
    ///    
    /// \param c [in] Scalar constant.
    /// \return Reference to the destination matrix after the mutliplication.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////        
    Matrix&	operator*= ( const double c ) 
    { 
      std::transform( mData, mData+mSize, mData, std::bind2nd(std::multiplies<double>(),c) ); 
      return(*this);  
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Scalar inplace division operator.
    ///    
    /// \param c [in] Scalar constant.
    /// \return Reference to the destination matrix after the division.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////        
	  Matrix& operator/= ( const double c ) 
    { 
      std::transform( mData, mData+mSize, mData, std::bind2nd(std::divides<double>(),c) ); 
      return(*this);
    }
           
    Matrix&	operator+= ( const SparseMatrix& m );	    
   	Matrix& operator-= ( const SparseMatrix& m );		            
    Matrix& operator*= ( const Matrix& m );            
    Matrix& operator*= ( const SparseMatrix& m );
  //@}    
    
  /// \name Accessors
  //@{
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief # of rows accessor.
    /// \return # of rows in the matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  size_t GetNumRows( void ) const 
    { return( mM ); }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief # of columns accessor.
    /// \return # of columns in the matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  size_t GetNumColumns( void ) const 
    { return( mN ); }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief %Matrix size assessor.
    /// \param M [out] # of rows in the matrix.
    /// \param N [out] # of columns in the matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void GetSize( size_t& M, size_t& N ) const 
    { M = mM; N = mN; }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief %Matrix element accessor by reference.
    /// \param i [in] Row index.
    /// \param j [in] Column index.
    /// \return Reference to the matrix element (i,j).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double& operator() ( const size_t i, const size_t j ) 
    { return( mData[ Index(i,j) ] ); }
	
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief %Matrix element accessor
    /// \param i [in] Row index.
    /// \param j [in] Column index.
    /// \return Value of the matrix element (i,j).
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double operator() (const size_t i, const size_t j ) const 
    { return( mData[ Index(i,j) ] ); }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Checks if the matrix is square.
    /// \return True if the # of rows equals the # of columns.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool IsSquare( void ) const 
    { return( mM == mN ); }
  //@}
       
  //@{ 
    /// \name Block assignment functions.
    void SetRow( const int i, const Vector& v );	
    void SetColumn( const int j, const Vector& v );
    void InsertBlock( const SparseMatrix& m, const int i, const int j );    
    void InsertBlock( const Matrix& m, const int i, const int j );
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Clears a row in the matrix.
    /// \param i [in] Index of the row to clear.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void ClearRow( const int i ) 
    { 
      for ( int j = 0; j < mN; ++j )
        mData[ Index(i,j) ] = 0.0;
    }
  //@}        
    
  //@{ 
    /// \name Column major array functions.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Builds a column major matrix array for use with other linear algebra libraries.    
    /// \warning The memory pointed to by the array is owned by the called and must be explicitly released.
    /// \warning Any memory pointed to by the array before the function call is not released by the function.
    /// \param A [out] Reference to a column major array that contains the matrix values.
    /// \return Pointer to the array.  It points to the same memory location as the parameter A after the function call.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////        
    double* GetFullStorageArray( double*& A ) const
    {
      A = new double[mSize];
      Assert( A != NULL, "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
      for ( int i = 0; i < mM; ++i )
        for ( int j = 0; j < mN; ++j )
          A[ j*mM + i ] = mData[ Index(i,j) ];
      
      return( A );
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Sets the matrix values equal to the values in a given column major matrix array.
    /// \warning It is assumed that the length of the array is the same as the number of elements in the matrix (MxN).
    /// \param A [in] Column major array that holds the matrix values.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////      
    void SetFullStorageArray( double* A )
    {
      Assert( A != NULL, "Passed a NULL Array", __FILE__, __FUNCTION__, __LINE__ );
       for ( int i = 0; i < mM; ++i )
        for ( int j = 0; j < mN; ++j )
          mData[ Index(i,j) ] = A[ j*mM + i ];
    }
  //@}
    
  //@{
    /// \name Miscellaneous
    void MakeTestMatrix( void );
  //@}

  //@{
    /// \name Stream operators
	  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Output stream insert operator.
    ///
    /// Outputs the matrix in tab delimted format.
    /// \param o [in, out] Output stream to insert the matrix values into.
    /// \param m [in] Matrix to insert into the stream.
    /// \return Reference to the output stream after the matrix values are inserted.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	   
    friend std::ostream& operator<< ( std::ostream& o, const Matrix& m ) 
    { 
      for ( int i = 0; i < m.mM; ++i ) 
      {
        std::copy( m.mData+m.Index(i,0), m.mData+m.Index(i+1,0), std::ostream_iterator<double>(o, "\t") );
        o << std::endl; 
      } 
      return(o); 
    }
  //@}

  //@{
    /// \name Elementry operators
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix-Matrix addition operator.
    /// \warning The sizes of the matrices must match.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    friend Matrix operator+ ( const Matrix& A, const Matrix& B ) 
    { 
      Matrix result(A); 
      return( result += B ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix-Matrix subtraction operator.
    /// \warning The sizes of the matrices must match.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator- ( const Matrix& A, const Matrix& B ) 
    { 
      Matrix result(A); 
      return( result -= B ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix-Matrix multiplication operator.
    /// \warning The number of columns of the first matrix must equal the number of rows of the second matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator* ( const Matrix& A, const Matrix& B ) 
    { 
      Matrix result(A); 
      return( result *= B ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix-Scalar multiplication operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator* ( const Matrix& A, const double c ) 
    { 
      Matrix result(A); 
      return( result *= c ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Scalar-Matrix multiplication operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator* ( const double c, const Matrix& A ) 
    { 
      Matrix result(A); 
      return( result *= c ); 
    } 
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix negation operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator- ( const Matrix& A ) 
    { 
      Matrix result(A); 
      return( result *= -1 ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix-SparseMatrix multiplication operator.
    /// \warning The dense matrix must be square and the same size as the sparse matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator* ( const Matrix& A, const SparseMatrix& B ) 
    { 
      Matrix result(A); 
      return( result *= B ); 
    }
        
    friend Matrix operator* ( const SparseMatrix& A, const Matrix& B );
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix-SparseMatrix addition operator
    /// \warning The dense matrix must be square and the same size as the sparse matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator+ ( const Matrix& A, const SparseMatrix& B ) 
    { 
      Matrix result(A); 
      return( result += B ); 
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief SparseMatrix-Matrix addition operator.
    /// \warning The dense matrix must be square and the same size as the sparse matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator+ ( const SparseMatrix& A, const Matrix& B ) 
    { 
      Matrix result(B);
      return( result += A ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix-SparseMatrix subtraction operator.
    /// \warning The dense matrix must be square and the same size as the sparse matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator- ( const Matrix& A, const SparseMatrix& B ) 
    { 
      Matrix result(A); 
      return( result -= B ); 
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief SparseMatrix-Matrix subtraction operator.
    /// \warning The dense matrix must be square and the same size as the sparse matrix.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator- ( const SparseMatrix& A, const Matrix& B ) 
    { 
      Matrix result(B); 
      return( -(result -= A) ); 
    }
        
    friend Vector operator* ( const Matrix& A, const Vector& x );
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Matrix-Scalar division operator.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    friend Matrix operator/ ( const Matrix& A, const double c )
    {
      Matrix result(A);
      return( result/c );
    }
  //@}
  
  //@{
    /// \name %Matrix functions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Computes the Frobenius norm of the matrix.
    /// 
    /// Computes the norm
    /// \f[ \| M \| = \sqrt{ \sum_{i,j} M_{i,j}^{2} } \f]
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    friend double frobnorm( const Matrix& M ) 
    { return( ::sqrt( std::inner_product(M.mData, M.mData+M.mSize, M.mData, 0.0) ) ); }
            
    friend double trace( const Matrix& M );                
    friend Matrix transpose( const Matrix& M );
    
    friend Matrix outerProduct( const Vector& u, const Vector& v );

    friend Matrix inv( const Matrix& M );
    friend Matrix exp( const Matrix& M );
    friend Matrix sqrt( const Matrix& M );
    friend Matrix invsqrt( const Matrix& M );
    friend Matrix cbrt( const Matrix& M );
    friend Matrix invcbrt( const Matrix& M );
    friend Matrix pow( const Matrix& M, const double c );
    friend Matrix ln( const Matrix& M );
    friend Matrix log10( const Matrix& M );
    friend Matrix cos( const Matrix& M );
    friend Matrix sin( const Matrix& M );
    friend void sincos( const Matrix& M, Matrix& C, Matrix& S );
    friend Matrix tan( const Matrix& M );
    friend Matrix acos( const Matrix& M );
    friend Matrix asin( const Matrix& M );
    friend Matrix atan( const Matrix& M );
    friend Matrix cosh( const Matrix& M );
    friend Matrix sinh( const Matrix& M );
    friend Matrix tanh( const Matrix& M );
    friend Matrix acosh( const Matrix& M );
    friend Matrix asinh( const Matrix& M );
    friend Matrix atanh( const Matrix& M );
    friend Matrix erf( const Matrix& M );
    friend Matrix erfc( const Matrix& M );

  //@}
    
  private:
    void AllocateMemory( void ) 
    {
      //Assert( mM >= 0, "Zero or Negative # of Rows", __FILE__, __FUNCTION__, __LINE__ ); 
      //Assert( mN >= 0, "Zero or Negative # of Columns", __FILE__, __FUNCTION__, __LINE__ );
      Assert( mData == NULL, "Data Array Already Allocated", __FILE__, __FUNCTION__, __LINE__ );           
      
      if ( mSize != 0 )  
      {
        mData = new double[mSize];
        Assert( mData != NULL, "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
      }       
    }
    
	  void ReleaseMemory( void ) 
    { SafeDeleteArray(mData); }
    
	  void StraightCopy( const Matrix& m )     
    { std::copy( m.mData, m.mData+m.mSize, mData ); }
    
	  void AllocateCopy( const Matrix& m ) 
    { 
      ReleaseMemory(); 
      AllocateMemory(); 
      StraightCopy(m); 
    }
    
    void StraightCopy( const SparseMatrix& m );
    
    void AllocateCopy( const SparseMatrix& m ) 
    { 
      ReleaseMemory(); 
      AllocateMemory();
      StraightCopy(m); 
    }
    
    bool SizesMatch( const size_t M, const size_t N ) const 
    { return( (mM == M) && (mN == N) ); }
    
	  size_t Index( const size_t i, const size_t j) const 
    { return( i*mN + j ); }
	
    size_t mM;
    size_t mN;
    size_t	mSize;  
    double* mData;
  };   
}

#endif

