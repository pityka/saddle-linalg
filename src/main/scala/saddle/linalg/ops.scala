package org.saddle.linalg

trait AxB
trait AtxB
trait AxBt
trait AtxBt

trait aAxBpbC
trait aAtxBpbC
trait aAxBtpbC
trait aAtxBtpbC

trait AxAt
trait AtxA
trait aAxAtpbC
trait aAtxApbC

trait InvertWithLU
trait InvertPDCholesky

trait GeneralSVD

// from lapack:
// inverse of A'A from A (QR) (http://scicomp.stackexchange.com/questions/3188/dealing-with-the-inverse-of-a-positive-definite-symmetric-covariance-matrix)
// svd https://en.wikipedia.org/wiki/Singular_value_decomposition#Calculating_the_SVD
// http://www.netlib.org/lapack/lapack-3.1.1/html/dgesdd.f.html
// eigen (symmetric + nonsymmetric)
