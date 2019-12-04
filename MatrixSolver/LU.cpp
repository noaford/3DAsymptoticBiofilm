#include "LU.hpp"

#include "LinearAlgebraSettings.hpp"

#ifdef __INTEL_COMPILER
  #include <omp.h>
#endif

// External Intel Math Kernel Library Direct Sparse Matrix Solver
#if 0
#define PARDISO pardiso_

extern "C" int PARDISO( void* , int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, int *, double *, double *, int *);
#else
#include "mkl_types.h"
#include "mkl_pardiso.h"
#define PARDISO pardiso
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Constructs a &LU factorization of a sparse matrix.
/// \param A [in] Sparse matrix to factorize.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::LU::LU( const SparseMatrix& A ) : mA(NULL), mIa(NULL), mJa(NULL)
{   
  int phase;  // Current Solver Phase
           
  double ddum;  // Double Dummy Variable
  int idum;     // Integer Dummy Variable
  
  memset( mIparm, 0, sizeof(mIparm) );    // Zero Out Parameter Array
  
  mIparm[ 0] = 1;                                                 // Do Not Use Solver Default
  mIparm[ 1] = 2;                                                 // Use Fill-in Reordering from METIS
  
  #ifndef __GNUC__
    mIparm[ 2] = omp_get_max_threads();                           // Number of Processors from OMP_NUM_THREADS
  #else
    mIparm[ 2] = 1;                                               // Single threaded for GCC 
  #endif
  
  mIparm[ 3] = 0;                                                 // Do Not Use Iterative-Direct Algorithm
  mIparm[ 4] = 0;                                                 // Do Not Use User Fill-in Reducing Permutation
  mIparm[ 5] = 0;                                                 // Write Solution into x
  mIparm[ 6] = 0;                                                 // Not In Use
  mIparm[ 7] = 0;                                                 // Maximum Number of Iterative Refinement Steps
  mIparm[ 8] = 0;                                                 // Not In Use
  mIparm[ 9] = LinearAlgebra::Settings::DirectSolverPerturb;      // Amount to Perturb the Pivot Elements
  mIparm[10] = 1;                                                 // Use Nonsymmetric Permutation and Scaling MPS
  mIparm[11] = 0;                                                 // Not In Use
  mIparm[12] = 0;                                                 // Not In Use
  mIparm[13] = -1;                                                // Output the # of Perturbed Pivots
  mIparm[14] = -1;                                                // Output the Peak Memory Usage in KB for Reordering
  mIparm[15] = -1;                                                // Output the Memory Required for Factorization
  mIparm[16] = -1;                                                // Output the Total Memory Used for Solving
  mIparm[17] = -1;                                                // Output the # of Nonzeros in the Factor LU
  mIparm[18] = -1;                                                // Output the # of Mflops for LU Factorization
  mIparm[19] = -1;                                                // Output the # of CG Iterations
  
  int mtype  = 11;  // Real Unsymmetic Matrix
  int nrhs   = 1;   // # of Right Hand Sides
  int maxfct = 1;   // Maximum # of Numerical Factorizations
  int mnum   = 1;   // Use the First Matrix
  int msglvl = 0;   // Do Not Print Statistical Information in a File
  int error  = 0;   // Initialize Error Flag
  
  memset( mPt, 0, sizeof(void*)*64 );  // Zero Out the Internal Solver Pointer  
        
  CreateMatrix( A ); 
    
  phase = 12;  // Reordering and Symbolic Factorization Phase
    MKL_INT mkl_mN = MKL_INT(mN);
  PARDISO( mPt, &maxfct, &mnum, &mtype, &phase, &mkl_mN, mA, mIa, mJa, &idum, &nrhs, mIparm, &msglvl, &ddum, &ddum, &error );
      
  if ( error != 0 )
  {
    ShowError( error );
    Assert( false, "Solver Error", __FILE__, __FUNCTION__, __LINE__ );
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Destructor.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
LinearAlgebra::LU::~LU( void )
{
  SafeDeleteArray( mA );
  SafeDeleteArray( mIa );
  SafeDeleteArray( mJa );
  
  int phase = -1;  // Release Internal Memory
  int msglvl = 0;
  int error = 0;
  int idum;
  double ddum;
  PARDISO( mPt, &idum, &idum, &idum, &phase, &idum, &ddum, &idum, &idum, &idum, &idum, mIparm, &msglvl, &ddum, &ddum, &error );
  
  if ( error != 0 )
    ShowError( error );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Solves a linear system \f$ L U \vec{x} = \vec{b} \f$ for \f$ \vec{x} \f$.
/// \warning The length of b must be the same size as the dimension of the matrix.
/// \warning The length of x must be the same size as the dimension of the matrix.
/// \param b [in] RHS vector.
/// \param x [in, out] Solution vector. 
/// \return True if the solve was successful.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool LinearAlgebra::LU::Solve( const Vector& b, Vector& x ) const 
{
  Assert( b.GetSize() == mN, "Matrix and Vector Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  Assert( x.GetSize() == mN, "Matrix and Vector Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
  
  int mtype  = 11;  // Real Unsymmetic Matrix
  int nrhs   = 1;   // # of Right Hand Sides
  int maxfct = 1;   // Maximum # of Numerical Factorizations
  int mnum   = 1;   // Use the First Matrix
  int msglvl = 0;   // Do Not Print Statistical Information in a File
  int error  = 0;   // Initialize Error Flag
  int idum;     // Integer Dummy Variable
  
  double* pB = NULL;
  double* pX = NULL;
      
  CreateVectors( b, pB, pX );
  
  int phase = 33;  // Back Substitution and Iterative Refinement 
    MKL_INT mkl_mN = MKL_INT(mN);
  PARDISO( const_cast<void**>(mPt), &maxfct, &mnum, &mtype, &phase, &mkl_mN, mA, mIa, mJa, &idum, &nrhs, const_cast<int*>(mIparm), &msglvl, pB, pX, &error );   
  
  if ( error != 0 )
  {
    ShowError( error );
    return( false );
  }
  
  x.SetFullStorageArray( pX );
  
  SafeDeleteArray(pB);
  SafeDeleteArray(pX);
  
  return( true );
}

void LinearAlgebra::LU::CreateMatrix( const SparseMatrix& A )
{ 
  A.GetCompressedFormatArrays( mA, mIa, mJa ); 
  mN = A.GetSize();
}

void LinearAlgebra::LU::CreateVectors( const Vector& b, double*& pB, double*& pX ) const
{  
  pX = new double[mN];
  Assert( pX != NULL, "Unable to Allocate Memory", __FILE__, __FUNCTION__, __LINE__ );
  
  b.GetFullStorageArray( pB );
}

void LinearAlgebra::LU::ShowError( const int error ) const
{
  std::cerr << "LU Decomposition Error - ";
    
  switch ( error )
    {
    case (-1):
      std::cerr << "Inconsistent Input ";
      break;
        
    case (-2):  
      std::cerr << "Not Enough Memory";
      break;
        
    case (-3):
      std::cerr << "Reordering Problem";
      break;
        
    case (-4):
      std::cerr << "Zero Pivot, Numerical Factorization Problem";
      break;
        
    case (-5):  
      std::cerr << "Unclassified (Internal) Error";
      break;
        
    case (-6):  
      std::cerr << "Perordering Failed";
      break;
        
    case (-7):
      std::cerr << "Diagonal Matrix Problem";
      break;
        
    default:
      std::cerr << "Unknown Error - " << error;
      break;
    };
    
  std::cerr << std::endl;
}
