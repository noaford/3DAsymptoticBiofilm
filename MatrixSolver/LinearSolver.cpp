#include "LinearSolver.hpp"

#include "SparseMatrix.hpp"
#include "LinearAlgebra.hpp"
#include "LinearAlgebraSettings.hpp"
#include "LongTimer.hpp"
#include "Vector.hpp"

#if defined(MKL_ILP64)
#define MKL_INT long long
#else
#define MKL_INT long
#endif

#ifdef __INTEL_COMPILER
#include <omp.h>
#endif

extern "C" void PARDISO( void * , MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, double *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, double *, double *, MKL_INT *);


std::string LinearAlgebra::DirectSolver::DirectSolverHeader("\t** Direct System Solver - ");

bool LinearAlgebra::DirectSolver::Solve( const SparseMatrix& matA, const SparseVector& vB, SparseVector& vX, bool Output )
{
    Assert( matA.GetSize() == vB.GetSize(), "Matrix and Right Hand Side Vector are Not the Same Size", __FILE__, __FUNCTION__, __LINE__ );
    Assert( matA.GetSize() == vX.GetSize(), "Matrix and Solution Vector are Not the Same Size", __FILE__, __FUNCTION__, __LINE__ );
    
    LongTimer timer;
    
    MKL_INT mtype = 11;  // Real Unsymmetic Matrix
    MKL_INT nrhs  = 1;   // # of Right Hand Sides
    
    void* pt[64];   // Internal Solver Memory Pointer
    
    MKL_INT iparm[64];  // Pardiso Control Parameters
    
    MKL_INT phase;  // Current Solver Phase
    
    double ddum;  // Double Dummy Variable
    MKL_INT idum;     // Integer Dummy Variable
    
    for(int i=0;i<64;i++)
        iparm[i] = 0;
    
    iparm[ 0] = 1;                                                 // Do Not Use Solver Default
    iparm[ 1] = 2;                                                 // Use Fill-in Reordering from METIS
    
#ifdef __INTEL_COMPILER
    iparm[ 2] = omp_get_max_threads();                           // Number of Processors from OMP_NUM_THREADS
#else
    iparm[ 2] = 1;                                               // Single threaded if not using Intel C++ compiler
#endif
    
    iparm[ 3] = 0;                                                 // Do Not Use Iterative-Direct Algorithm
    iparm[ 4] = 0;                                                 // Do Not Use User Fill-in Reducing Permutation
    iparm[ 5] = 0;                                                 // Write Solution into x
    iparm[ 6] = 0;                                                 // Not In Use
    iparm[ 7] = 0;                                                 // Maximum Number of Iterative Refinement Steps
    iparm[ 8] = 0;                                                 // Not In Use
    iparm[ 9] = LinearAlgebra::Settings::DirectSolverPerturb;      // Amount to Perturb the Pivot Elements
    iparm[10] = 1;                                                 // Use Nonsymmetric Permutation and Scaling MPS
    iparm[11] = 0;                                                 // Not In Use
    iparm[12] = 0;                                                 // Not In Use
    iparm[13] = -1;                                                // Output the # of Perturbed Pivots
    iparm[14] = -1;                                                // Output the Peak Memory Usage in KB for Reordering
    iparm[15] = -1;                                                // Output the Memory Required for Factorization
    iparm[16] = -1;                                                // Output the Total Memory Used for Solving
    iparm[17] = -1;                                                // Output the # of Nonzeros in the Factor LU
    iparm[18] = -1;                                                // Output the # of Mflops for LU Factorization
    iparm[19] = -1;                                                // Output the # of CG Iterations
    
    
    
    MKL_INT maxfct = 1;                                                // Maximum # of Numerical Factorizations
    MKL_INT mnum   = 1;                                                // Use the First Matrix
    MKL_INT msglvl = LinearAlgebra::Settings::DirectSolverFileOutput;  // Print Statistical Information in a File
    MKL_INT error  = 0;                                                // Initialize Error Flag
    
    memset( pt, 0, sizeof(void*)*64 );  // Zero Out the Internal Solver Pointer
    
    
    
    double* A = NULL;  // Sparse Matrix Values
    double* b = NULL;  // RHS Vector
    double* x = NULL;  // Solution Vector
    
    const double nSize = matA.GetSize();
    const double nNonZero = matA.GetNumNonZeroElements();
    if ( Output )
    {
        std::cout << DirectSolverHeader << "System Size = " << nSize << std::endl;
        std::cout << DirectSolverHeader << "# of Nonzeros = " << nNonZero << std::endl;
        std::cout << DirectSolverHeader << "% Sparse = " << 100.0*(nNonZero/nSize/nSize) << "%" << std::endl;
        
        std::cout << DirectSolverHeader << "Assembling . . . ";
        std::cout.flush();
    }
    size_t intn;
    int* intia = NULL;
    int* intja = NULL;
    
    AssembleArrays( matA, vB, intn, intia, intja, A, b, x );  // Assemble Arrays for the Solver from Sparse Matrix and Vector
    
    //Convert int into MKL_INT
    MKL_INT n = static_cast<MKL_INT>(intn);
    
    size_t lengtharrayia = intn+1;
    MKL_INT* ia = new MKL_INT[lengtharrayia];
    
    for(int i=0; i<lengtharrayia; i++)
    {
        ia[i] = static_cast<MKL_INT>(intia[i]);
    }
    
    int lengtharrayja = nNonZero;
    MKL_INT* ja = new MKL_INT[lengtharrayja];
    for(int i=0; i<lengtharrayja; i++)
    {
        ja[i] = static_cast<MKL_INT>(intja[i]);
    }
    if ( Output )
    {
        std::cout << "Completed. " << std::endl;
        std::cout << DirectSolverHeader << "Reordering . . . ";
        std::cout.flush();
    }
    phase = 11;  // Reordering and Symbolic Factorization Phase
    
    PARDISO( pt, &maxfct, &mnum, &mtype, &phase, &n, A, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error );
    
    if ( error != 0 )
    {
        ShowDirectError( static_cast<int>(error) );
        return( false );
    }
    
    if ( Output )
    {
        std::cout << "Completed. " << std::endl;
        std::cout << DirectSolverHeader << "# of Nonzeros in LU Factors = " << iparm[17] << std::endl;
        std::cout << DirectSolverHeader << "# of Factorization MFLOPS = " << iparm[18] << std::endl;
        //std::cout << DirectSolverHeader << "Peak Memory Usage for Reordering = " << std::setprecision(4) << iparm[14]/1024.0 << " MB " << std::endl;
        //std::cout << DirectSolverHeader << "Memory Needed for Numerical Factorization = " << std::setprecision(4) << iparm[15]/1024.0 << " MB " << std::endl;
        std::cout << DirectSolverHeader << "Factoring . . . ";
        std::cout.flush();
    }
    
    phase = 22;  // Numerical Factorization Phase
    PARDISO( pt, &maxfct, &mnum, &mtype, &phase, &n, A, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error );
    
    if ( error != 0 )
    {
        ShowDirectError( static_cast<int>(error) );
        return( false );
    }
    
    if ( Output )
    {
        std::cout << "Completed. " << std::endl;
        //std::cout << DirectSolverHeader << "Total Memory Used = " << std::setprecision(4) << iparm[16]/1024.0 << " MB" << std::endl;
        std::cout << DirectSolverHeader << "# of Perturbed Pivots = " << iparm[13] << std::endl;
        // std::cout << DirectSolverHeader << "# of CGS Iterations = ";
        // if (iparm[19] >= 0 )
        // std::cout << iparm[19] << std::endl;
        // else
        // {
        // std::cout << -iparm[19]/10 << std::endl;
        // std::cerr << DirectSolverHeader << "ERROR (CGS) - ";
        // switch( -iparm[19] % 10 )
        // {
        // case 1:
        // std::cerr << "Too Large Fluctuations of the Residuum" << std::endl;
        // break;
        // case 2:
        // std::cerr << "|| Dx_(max_it_cgs/2) || Too Large (Slow Convergence)" << std::endl;
        // break;
        // case 3:
        // std::cerr << "Stopping Criterion Not Reached at Maximum Iterations" << std::endl;
        // break;
        // case 4:
        // std::cerr << "Perturbed Pivots Caused Iterative Refinement" << std::endl;
        // break;
        // case 5:
        // std::cerr << "Factorization Is Too Fast for this Matrix" << std::endl;
        // break;
        // default:
        // std::cerr << "Unknown Error (" << (-iparm[19] % 10) << ")" << std::endl;
        // break;
        // }
        // }
        std::cout << DirectSolverHeader << "Backsolving . . . " << std::flush;
    }
    
    phase = 33;  // Back Substitution and Iterative Refinement
    PARDISO( pt, &maxfct, &mnum, &mtype, &phase, &n, A, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error );
    
    if ( Output )
        std::cout << "Completed. " << std::endl;
    
    if ( error != 0 )
    {
        ShowDirectError( static_cast<int>(error) );
        return( false );
    }
    
    if ( Output )
    {
        std::cout << DirectSolverHeader << "Assembling Solution . . . ";
        std::cout.flush();
    }
    
    AssembleVector( vX, x );
    
    if ( Output )
        std::cout << "Completed." << std::endl;
    
    phase = -1;  // Release Internal Memory
    
    if ( Output )
    {
        std::cout << DirectSolverHeader << "Releasing Internal Memory . . . ";
        std::cout.flush();
    }
    
    PARDISO( pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error );
    
    if ( Output )
    {
        std::cout << "Completed." << std::endl;
        std::cout << DirectSolverHeader << "Solve Completed in " << timer.HrsMinSec() << " seconds. " << std::endl;
    }
    
    // Release Arrays
    SafeDeleteArray( ia );
    SafeDeleteArray( ja );
    SafeDeleteArray( A );
    SafeDeleteArray( b );
    SafeDeleteArray( x );
    
    return( true );
}

void LinearAlgebra::DirectSolver::AssembleArrays( const SparseMatrix& matA, const SparseVector& vB, size_t& n, int*& ia, int*& ja, double*& A, double*& b, double*& x )
{
    // Make Sure the Arrays are Not Already Allocated and Release Them if They Are
    SafeDeleteArray( ia );
    SafeDeleteArray( ja );
    SafeDeleteArray( A );
    SafeDeleteArray( b );
    SafeDeleteArray( x );
    
    // Assert the System Matrix and Vector are the Same Size
    Assert( matA.GetSize() == vB.GetSize(), "Matrix and Vector Sizes Do Not Match", __FILE__, __FUNCTION__, __LINE__ );
    
    // Set the System Size;
    n = matA.GetSize();
    
    // Allocate Solution Array
    x  = new double[n];
    Assert( x != NULL, "Unable to Allocate Solution Array", __FILE__, __FUNCTION__, __LINE__ );
    
    // Copy the Vector vB into RHS Array
    vB.GetFullStorageArray( b );
    
    // Fill the Array of Row Indices, Column Indices, and Matrix Values
    matA.GetCompressedFormatArrays( A, ia, ja );
}

void LinearAlgebra::DirectSolver::AssembleVector( SparseVector& vX, double* x )
{ vX.SetFullStorageArray( x ); }

void LinearAlgebra::DirectSolver::ShowDirectError( const int error )
{
    std::cerr << DirectSolverHeader << " ERROR - ";
    
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
            std::cerr << "Preordering Failed";
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
