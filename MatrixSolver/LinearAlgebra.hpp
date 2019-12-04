#ifndef _LINEARALGEBRA_HPP_10142004_
#define _LINEARALGEBRA_HPP_10142004_

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Forward declarations for linear algebra classes 
/// \ingroup linalg
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace LinearAlgebra
{
  class Matrix;
  class Vector;
  class SparseMatrix;
  class SparseMatrixBlock;
  class DiagonalMatrix;
  class Tensor3;
  
  typedef Vector SparseVector;

  DiagonalMatrix inverse(const DiagonalMatrix& A);
  double trace(const DiagonalMatrix& A);
  double frobnorm(const DiagonalMatrix& A);
  
  Matrix outerProduct(const Vector& u, const Vector& v);
  Matrix exp(const Matrix& M);
  Matrix sqrt(const Matrix& M);
  Matrix cbrt(const Matrix& M);
  Matrix pow(const Matrix& M, const double c);
  Matrix ln(const Matrix& M);
  Matrix log10(const Matrix& M);
  Matrix cos(const Matrix& M);
  Matrix sin(const Matrix& M);
  void sincos(const Matrix& M, Matrix& C, Matrix& S);
  Matrix tan(const Matrix& M);
  Matrix acos(const Matrix& M);
  Matrix asin(const Matrix& M);
  Matrix atan(const Matrix& M);
  Matrix cosh(const Matrix& M);
  Matrix sinh(const Matrix& M);
  Matrix tanh(const Matrix& M);
  Matrix acosh(const Matrix& M);
  Matrix asinh(const Matrix& M);
  Matrix atanh(const Matrix& M);
  Matrix erf(const Matrix& M);
  Matrix erfc(const Matrix& M);

  Vector trans_mult(const SparseMatrix& M, const Vector& v);
  double infnorm(const SparseMatrix& M);
  double onenorm(const SparseMatrix& M);
  double cond(const SparseMatrix& M);
  double RowMult(const SparseMatrix& M, const Vector& v, int i);
  double maxEV(const SparseMatrix& M, Vector& v, double a, int i);
  double minEV(const SparseMatrix& M, Vector& v, double a, int i);
}

#endif
