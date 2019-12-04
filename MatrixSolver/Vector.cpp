#include "Vector.hpp"

#include "SparseMatrix.hpp"

#include "GNUPlot.hpp"

LinearAlgebra::Vector& LinearAlgebra::Vector::operator*= ( const LinearAlgebra::SparseMatrix& A )
{
  (*this) = A*(*this);

  return(*this);
}

void LinearAlgebra::Plot( const LinearAlgebra::Vector& a, const std::string& Title )
{
  //const int nSize = a.GetSize();
  
  //double  pX[nSize];   
  
  //for ( int i = 0; i <= nSize; ++i )
  //  pX[i] = a(i); 
  
  //GNUPlot plot;
     
  //plot.SetStyle( "lines" );
  //plot.Plot( pX, nSize, Title.c_str() );
      
  //return;
}
