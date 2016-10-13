package org.saddle.linalg

import org.saddle._
import scala.{specialized => spec}
import annotation.implicitNotFound

class Pimp(val self: Mat[Double]) extends LinalgOps

@implicitNotFound(msg = "${O} not found")
trait MatBinOp[O, Res] {
  def apply(a: Mat[Double], b: Mat[Double]): Res
}

@implicitNotFound(msg = "${O} not found")
trait MatGemmOp[O, Res] {
  def apply(a: Mat[Double],
            b: Mat[Double],
            c: Mat[Double],
            alpha: Double,
            beta: Double): Res
}

@implicitNotFound(msg = "${O} not found")
trait MatGemmSelfOp[O, Res] {
  def apply(a: Mat[Double], c: Mat[Double], alpha: Double, beta: Double): Res
}

@implicitNotFound(msg = "${O} not found")
trait MatUnaryOp[O, Res] {
  def apply(a: Mat[Double]): Res
}

@implicitNotFound(msg = "${O} not found")
trait MatUnaryOp1Scalar[O, T, Res] {
  def apply(a: Mat[Double], s: T): Res
}

trait LinalgOps {
  val self: Mat[Double]
  type B = Mat[Double]

  def linalg = this

  def invert(implicit op: MatUnaryOp[InvertWithLU, B]): B = op(self)

  def invertPD(implicit op: MatUnaryOp[InvertPDCholesky, B]): B = op(self)

  /**
    * Simple DGEMM
    */
  /* A x B */
  def mm(other: B)(implicit op: MatBinOp[AxB, B]): B =
    op(self, other)

  /* t(A) x B */
  def tmm(other: B)(implicit op: MatBinOp[AtxB, B]): B =
    op(self, other)

  /* A x t(B) */
  def mmt(other: B)(implicit op: MatBinOp[AxBt, B]): B =
    op(self, other)

  /* t(A) x t(B) */
  def tmmt(other: B)(implicit op: MatBinOp[AtxBt, B]): B =
    op(self, other)

  /**
    * Full DGEMM
    */
  /* alhpa A x B + beta * C */
  def mmc(other: B, c: B, alpha: Double = 1.0, beta: Double = 1.0)(
      implicit op: MatGemmOp[aAxBpbC, B]): B =
    op(self, other, c, alpha, beta)

  /* alhpa t(A) x B + beta * C */
  def tmmc(other: B, c: B, alpha: Double = 1.0, beta: Double = 1.0)(
      implicit op: MatGemmOp[aAtxBpbC, B]): B =
    op(self, other, c, alpha, beta)

  /* alhpa A x t(B) + beta * C */
  def mmtc(other: B, c: B, alpha: Double = 1.0, beta: Double = 1.0)(
      implicit op: MatGemmOp[aAxBtpbC, B]): B =
    op(self, other, c, alpha, beta)

  /* alhpa t(A) x t(B) + beta * C */
  def tmmtc(other: B, c: B, alpha: Double = 1.0, beta: Double = 1.0)(
      implicit op: MatGemmOp[aAtxBtpbC, B]): B =
    op(self, other, c, alpha, beta)

  def innerM(implicit op: MatUnaryOp[AtxA, Mat[Double]]): B = op(self)

  def innerMpC(alpha: Double, beta: Double, c: Mat[Double])(
      implicit op: MatGemmSelfOp[aAtxApbC, Mat[Double]]): B =
    op(self, c, alpha, beta)

  def outerM(implicit op: MatUnaryOp[AxAt, Mat[Double]]): B = op(self)

  def outerMpC(alpha: Double, beta: Double, c: Mat[Double])(
      implicit op: MatGemmSelfOp[aAxAtpbC, Mat[Double]]): B =
    op(self, c, alpha, beta)

  def svd(implicit op: MatUnaryOp[GeneralSVD, SVDResult]): SVDResult = op(self)

  def svd(max: Int)(
      implicit op: MatUnaryOp1Scalar[GeneralSVDTrunc, Int, SVDResult])
    : SVDResult = op(self, max)

  def trace(implicit op: MatUnaryOp[Trace, Double]): Double = op(self)

  def diag(implicit op: MatUnaryOp[Diag, Vec[Double]]): Vec[Double] = op(self)

  def isPositiveDefinite(implicit op: MatUnaryOp[TestPD, Boolean]): Boolean =
    op(self)

  def eigNonSymm(
      implicit op: MatUnaryOp[EigNS, EigenDecompositionNonSymmetric])
    : EigenDecompositionNonSymmetric =
    op(self)

  def eigSymm(implicit op: MatUnaryOp[EigS, EigenDecompositionSymmetric])
    : EigenDecompositionSymmetric =
    op(self)

  def eigSymm(i: Int)(
      implicit op: MatUnaryOp1Scalar[EigSTrunc,
                                     Int,
                                     EigenDecompositionSymmetric])
    : EigenDecompositionSymmetric =
    op(self, i)

  def solve(other: B)(implicit op: MatBinOp[GeneralSolve, B]): B =
    op(self, other)

}
