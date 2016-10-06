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

trait LinalgOps {
  val self: Mat[Double]
  type B = Mat[Double]

  def linalg = this

  /**
    * Simple DGEMM
    */
  /* A x B */
  def mm[Res](other: B)(implicit op: MatBinOp[AxB, Res]): Res =
    op(self, other)

  /* t(A) x B */
  def tmm[Res](other: B)(implicit op: MatBinOp[AtxB, Res]): Res =
    op(self, other)

  /* A x t(B) */
  def mmt[Res](other: B)(implicit op: MatBinOp[AxBt, Res]): Res =
    op(self, other)

  /* t(A) x t(B) */
  def tmmt[Res](other: B)(implicit op: MatBinOp[AtxBt, Res]): Res =
    op(self, other)

  /**
    * Full DGEMM
    */
  /* alhpa A x B + beta * C */
  def mmc[Res](other: B, c: B, alpha: Double = 1.0, beta: Double = 1.0)(
      implicit op: MatGemmOp[aAxBpbC, Res]): Res =
    op(self, other, c, alpha, beta)

  /* alhpa t(A) x B + beta * C */
  def tmmc[Res](other: B, c: B, alpha: Double = 1.0, beta: Double = 1.0)(
      implicit op: MatGemmOp[aAtxBpbC, Res]): Res =
    op(self, other, c, alpha, beta)

  /* alhpa A x t(B) + beta * C */
  def mmtc[Res](other: B, c: B, alpha: Double = 1.0, beta: Double = 1.0)(
      implicit op: MatGemmOp[aAxBtpbC, Res]): Res =
    op(self, other, c, alpha, beta)

  /* alhpa t(A) x t(B) + beta * C */
  def tmmtc[Res](other: B, c: B, alpha: Double = 1.0, beta: Double = 1.0)(
      implicit op: MatGemmOp[aAtxBtpbC, Res]): Res =
    op(self, other, c, alpha, beta)

}
