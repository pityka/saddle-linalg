package org.saddle.linalg

import org.saddle._

case class DPotrfException(i: Int)
    extends java.lang.Exception(s"""|dpotrf error, info=$i
              |*  INFO    (output) INTEGER
|*          = 0:  successful exit
|*          < 0:  if INFO = -i, the i-th argument had an illegal value
|*          > 0:  if INFO = i, the leading minor of order i is not
|*                positive definite, and the factorization could not be
|*                completed.""".stripMargin)

case class SVDResult(u: Mat[Double], sigma: Vec[Double], vt: Mat[Double])

trait OpImpl {

  lazy val BLAS = com.github.fommil.netlib.BLAS.getInstance

  lazy val LAPACK = com.github.fommil.netlib.LAPACK.getInstance

  implicit def svd = new MatUnaryOp[GeneralSVD, SVDResult] {
    def apply(m: Mat[Double]): SVDResult = {

      /** Lapack gives us the SVD of the transpose
        * t(a) = v t(s) t(u)
        *   a  = u s t(v)
        */
      val cop = m.contents.clone
      val s = Array.ofDim[Double](math.min(m.numRows, m.numCols))
      val u = Array.ofDim[Double](m.numCols * m.numCols)
      val vt = Array.ofDim[Double](m.numRows * m.numRows)
      val lworkQuery = Array.ofDim[Double](1)
      val info = new org.netlib.util.intW(0)

      // 1. Workspace query
      LAPACK.dgesvd(
        "A", // JOBU,
        "A", // JOBVT,
        m.numCols, // M,
        m.numRows, // N,
        cop, // A,
        m.numCols, // LDA,
        s, // S,
        u, // U,
        m.numCols, // LDU,
        vt, // VT,
        m.numRows, // LDVT,
        lworkQuery, // WORK,
        -1, // LWORK,
        info // INFO
      )

      val lwork = lworkQuery(0).toInt
      val work = Array.ofDim[Double](lwork)

      LAPACK.dgesvd(
        "A", // JOBU,
        "A", // JOBVT,
        m.numCols, // M,
        m.numRows, // N,
        cop, // A,
        m.numCols, // LDA,
        s, // S,
        u, // U,
        m.numCols, // LDU,
        vt, // VT,
        m.numRows, // LDVT,
        work, // WORK,
        lwork, // LWORK,
        info // INFO
      )

      if (lapackInfoMethod.get(info) == 0) {
        val ut: Mat[Double] = new mat.MatDouble(m.numRows, m.numRows, vt)
        val v: Mat[Double] = new mat.MatDouble(m.numCols, m.numCols, u)
        val sigma: Vec[Double] = Vec(s: _*)

        SVDResult(
          u = ut,
          sigma = sigma,
          vt = v
        )
      } else throw new RuntimeException("SVD Failed")

    }
  }

  implicit def invertGeneralLU = new MatUnaryOp[InvertWithLU, Mat[Double]] {
    def apply(m: Mat[Double]): Mat[Double] = {

      val marray = m.contents
      val array = marray.clone

      val ipiv = Array.ofDim[Int](math.max(1, math.min(m.numCols, m.numRows)))

      LAPACK.dgetrf(m.numCols,
                    m.numRows,
                    array,
                    m.numCols,
                    ipiv,
                    new org.netlib.util.intW(0))

      val lworkQuery = Array.ofDim[Double](1)

      LAPACK.dgetri(m.numCols,
                    array,
                    m.numCols,
                    ipiv,
                    lworkQuery,
                    -1,
                    new org.netlib.util.intW(0))

      val work = Array.ofDim[Double](lworkQuery(0).toInt + 1)
      LAPACK.dgetri(m.numCols,
                    array,
                    m.numCols,
                    ipiv,
                    work,
                    lworkQuery(0).toInt + 1,
                    new org.netlib.util.intW(0))

      new mat.MatDouble(m.numCols, m.numCols, array)

    }
  }

  private val lapackInfoMethod =
    java.lang.Class.forName("org.netlib.util.intW").getField("val")

  implicit def invertPD = new MatUnaryOp[InvertPDCholesky, Mat[Double]] {
    def apply(m: Mat[Double]): Mat[Double] = {

      val marray = m.contents
      val array = marray.clone
      val info = new org.netlib.util.intW(0)
      val info2 = new org.netlib.util.intW(0)

      LAPACK.dpotrf("L", m.numCols, array, m.numCols, info)
      LAPACK.dpotri("L", m.numCols, array, m.numCols, info2)

      if (lapackInfoMethod.get(info) == 0 && lapackInfoMethod
            .get(info2) == 0) {

        var i = 0
        var j = 0
        while (i < m.numCols) {
          while (j < i) {
            array(i * m.numCols + j) = array(j * m.numCols + i)
            j += 1
          }
          j = 0
          i += 1
        }

        new mat.MatDouble(m.numCols, m.numCols, array)

      } else if (lapackInfoMethod.get(info) != 0) {
        throw DPotrfException(lapackInfoMethod.get(info).asInstanceOf[Int])
      } else {
        throw new RuntimeException(
          "ERROR in dpotri info=" + lapackInfoMethod
            .get(info2) + """
                      |lapack says:
                      |      INFO    (output) INTEGER
                    |= 0:  successful exit
                    |< 0:  if INFO = -i, the i-th argument had an illegal
                    |value
                    |> 0:  if INFO = i, the (i,i) element of the factor U
                    |or L is zero, and the inverse could not be computed.""".stripMargin + ", matrix: " + m.toString)
      }

    }
  }

  implicit def mult1 =
    new MatBinOp[AxB, Mat[Double]] {
      def apply(a: Mat[Double], b: Mat[Double]): Mat[Double] = {
        assert(a.numCols == b.numRows,
               s"Incorrect dimensions ${a.numCols} ${b.numRows}")

        val result = Array.ofDim[Double](a.numRows * b.numCols)

        BLAS.dgemm("N",
                   "N",
                   b.numCols, // M
                   a.numRows, // N
                   b.numRows, // K
                   1.0, // alpha
                   b.contents, // op(a) data
                   b.numCols, // lda
                   a.contents, // op(b) data
                   b.numRows, //ldb
                   0.0, // c
                   result, // c data
                   b.numCols) // ldc
        new mat.MatDouble(a.numRows, b.numCols, result)
      }
    }

  implicit def mult1c =
    new MatGemmOp[aAxBpbC, Mat[Double]] {
      def apply(a: Mat[Double],
                b: Mat[Double],
                c: Mat[Double],
                alpha: Double,
                beta: Double): Mat[Double] = {
        assert(a.numCols == b.numRows,
               s"Incorrect dimensions ${a.numCols} ${b.numRows}")
        assert(c.numRows == a.numRows && c.numCols == b.numCols)

        val result = c.contents.clone

        BLAS.dgemm("N",
                   "N",
                   b.numCols, // M
                   a.numRows, // N
                   b.numRows, // K
                   alpha, // alpha
                   b.contents, // op(a) data
                   b.numCols, // lda
                   a.contents, // op(b) data
                   b.numRows, //ldb
                   beta, // c
                   result, // c data
                   b.numCols) // ldc
        new mat.MatDouble(a.numRows, b.numCols, result)
      }
    }

  implicit def mult2 =
    new MatBinOp[AtxB, Mat[Double]] {
      def apply(a: Mat[Double], b: Mat[Double]): Mat[Double] = {
        assert(a.numRows == b.numRows,
               s"Incorrect dimensions ${a.numRows} ${b.numRows}")

        val result = Array.ofDim[Double](a.numCols * b.numCols)

        BLAS.dgemm("N", // op a
                   "T", // op b
                   b.numCols, // M rows of op(a)
                   a.numCols, // N cols of op(b)
                   b.numRows, // K cols of op(a)
                   1.0, // alpha
                   b.contents, // op(a) data
                   b.numCols, // lda
                   a.contents, // op(b) data
                   a.numCols, //ldb
                   0.0, // c
                   result, // c data
                   b.numCols) // ldc
        new mat.MatDouble(a.numCols, b.numCols, result)
      }
    }

  implicit def mult2c =
    new MatGemmOp[aAtxBpbC, Mat[Double]] {
      def apply(a: Mat[Double],
                b: Mat[Double],
                c: Mat[Double],
                alpha: Double,
                beta: Double): Mat[Double] = {
        assert(a.numRows == b.numRows,
               s"Incorrect dimensions ${a.numRows} ${b.numRows}")
        assert(c.numRows == a.numCols && c.numCols == b.numCols)

        val result = c.contents.clone

        BLAS.dgemm("N", // op a
                   "T", // op b
                   b.numCols, // M rows of op(a)
                   a.numCols, // N cols of op(b)
                   b.numRows, // K cols of op(a)
                   alpha, // alpha
                   b.contents, // op(a) data
                   b.numCols, // lda
                   a.contents, // op(b) data
                   a.numCols, //ldb
                   beta, // c
                   result, // c data
                   b.numCols) // ldc
        new mat.MatDouble(a.numCols, b.numCols, result)
      }
    }

  implicit def mult2self =
    new MatUnaryOp[AtxA, Mat[Double]] {
      def apply(a: Mat[Double]): Mat[Double] = {

        val result = Array.ofDim[Double](a.numCols * a.numCols)

        BLAS.dgemm("N", // op a
                   "T", // op b
                   a.numCols, // M rows of op(a)
                   a.numCols, // N cols of op(b)
                   a.numRows, // K cols of op(a)
                   1.0, // alpha
                   a.contents, // op(a) data
                   a.numCols, // lda
                   a.contents, // op(b) data
                   a.numCols, //ldb
                   0.0, // c
                   result, // c data
                   a.numCols) // ldc
        new mat.MatDouble(a.numCols, a.numCols, result)
      }
    }

  implicit def mult2cself =
    new MatGemmSelfOp[aAtxApbC, Mat[Double]] {
      def apply(a: Mat[Double],
                c: Mat[Double],
                alpha: Double,
                beta: Double): Mat[Double] = {
        assert(c.numRows == a.numCols && c.numCols == a.numCols)

        val result = c.contents.clone

        BLAS.dgemm("N", // op a
                   "T", // op b
                   a.numCols, // M rows of op(a)
                   a.numCols, // N cols of op(b)
                   a.numRows, // K cols of op(a)
                   alpha, // alpha
                   a.contents, // op(a) data
                   a.numCols, // lda
                   a.contents, // op(b) data
                   a.numCols, //ldb
                   beta, // c
                   result, // c data
                   a.numCols) // ldc
        new mat.MatDouble(a.numCols, a.numCols, result)
      }
    }

  implicit def mult3 =
    new MatBinOp[AxBt, Mat[Double]] {
      def apply(a: Mat[Double], b: Mat[Double]): Mat[Double] = {
        assert(a.numCols == b.numCols,
               s"Incorrect dimensions ${a.numCols} ${b.numCols}")

        val result = Array.ofDim[Double](a.numRows * b.numRows)

        BLAS.dgemm("T",
                   "N",
                   b.numRows, // M
                   a.numRows, // N
                   b.numCols, // K
                   1.0, // alpha
                   b.contents, // op(a) data
                   b.numCols, // lda
                   a.contents, // op(b) data
                   b.numCols, //ldb
                   0.0, // c
                   result, // c data
                   b.numRows) // ldc
        new mat.MatDouble(a.numRows, b.numRows, result)
      }
    }

  implicit def mult3self =
    new MatUnaryOp[AxAt, Mat[Double]] {
      def apply(a: Mat[Double]): Mat[Double] = {

        val result = Array.ofDim[Double](a.numRows * a.numRows)

        BLAS.dgemm("T",
                   "N",
                   a.numRows, // M
                   a.numRows, // N
                   a.numCols, // K
                   1.0, // alpha
                   a.contents, // op(a) data
                   a.numCols, // lda
                   a.contents, // op(b) data
                   a.numCols, //ldb
                   0.0, // c
                   result, // c data
                   a.numRows) // ldc
        new mat.MatDouble(a.numRows, a.numRows, result)
      }
    }

  implicit def mult3c =
    new MatGemmOp[aAxBtpbC, Mat[Double]] {
      def apply(a: Mat[Double],
                b: Mat[Double],
                c: Mat[Double],
                alpha: Double,
                beta: Double): Mat[Double] = {
        assert(a.numCols == b.numCols,
               s"Incorrect dimensions ${a.numCols} ${b.numCols}")
        assert(c.numRows == a.numRows && c.numCols == b.numRows)

        val result = c.contents.clone

        BLAS.dgemm("T",
                   "N",
                   b.numRows, // M
                   a.numRows, // N
                   b.numCols, // K
                   alpha, // alpha
                   b.contents, // op(a) data
                   b.numCols, // lda
                   a.contents, // op(b) data
                   b.numCols, //ldb
                   beta, // c
                   result, // c data
                   b.numRows) // ldc
        new mat.MatDouble(a.numRows, b.numRows, result)
      }
    }

  implicit def mult3selfplus =
    new MatGemmSelfOp[aAxAtpbC, Mat[Double]] {
      def apply(a: Mat[Double],
                c: Mat[Double],
                alpha: Double,
                beta: Double): Mat[Double] = {
        assert(c.numRows == a.numRows && c.numCols == a.numRows)

        val result = c.contents.clone

        BLAS.dgemm("T",
                   "N",
                   a.numRows, // M
                   a.numRows, // N
                   a.numCols, // K
                   alpha, // alpha
                   a.contents, // op(a) data
                   a.numCols, // lda
                   a.contents, // op(b) data
                   a.numCols, //ldb
                   beta, // c
                   result, // c data
                   a.numRows) // ldc
        new mat.MatDouble(a.numRows, a.numRows, result)
      }
    }

  implicit def mult4 =
    new MatBinOp[AtxBt, Mat[Double]] {
      def apply(a: Mat[Double], b: Mat[Double]): Mat[Double] = {
        assert(a.numRows == b.numCols,
               s"Incorrect dimensions ${a.numRows} ${b.numCols}")

        val result = Array.ofDim[Double](a.numCols * b.numRows)

        BLAS.dgemm("T", // op a
                   "T", // op b
                   b.numRows, // M rows of op(a)
                   a.numCols, // N cols of op(b)
                   b.numCols, // K cols of op(a)
                   1.0, // alpha
                   b.contents, // op(a) data
                   b.numCols, // lda
                   a.contents, // op(b) data
                   a.numCols, //ldb
                   0.0, // c
                   result, // c data
                   b.numRows) // ldc
        new mat.MatDouble(a.numCols, b.numRows, result)
      }
    }

  implicit def mult4c =
    new MatGemmOp[aAtxBtpbC, Mat[Double]] {
      def apply(a: Mat[Double],
                b: Mat[Double],
                c: Mat[Double],
                alpha: Double,
                beta: Double): Mat[Double] = {
        assert(a.numRows == b.numCols,
               s"Incorrect dimensions ${a.numRows} ${b.numCols}")
        assert(c.numRows == a.numCols && c.numCols == b.numRows)
        val result = c.contents.clone

        BLAS.dgemm("T", // op a
                   "T", // op b
                   b.numRows, // M rows of op(a)
                   a.numCols, // N cols of op(b)
                   b.numCols, // K cols of op(a)
                   alpha, // alpha
                   b.contents, // op(a) data
                   b.numCols, // lda
                   a.contents, // op(b) data
                   a.numCols, //ldb
                   beta, // c
                   result, // c data
                   b.numRows) // ldc
        new mat.MatDouble(a.numCols, b.numRows, result)
      }
    }

}
