package org.saddle.linalg

import org.saddle._

class DPotrfException(i: Int) extends Exception(s"""|dpotrf error, info=$i
              |*  INFO    (output) INTEGER
|*          = 0:  successful exit
|*          < 0:  if INFO = -i, the i-th argument had an illegal value
|*          > 0:  if INFO = i, the leading minor of order i is not
|*                positive definite, and the factorization could not be
|*                completed.""".stripMargin)

case class SVDResult(u: Mat[Double], sigma: Vec[Double], vt: Mat[Double])

case class EigenDecompositionNonSymmetric(q: Mat[Double],
                                          lambdaReal: Vec[Double],
                                          lambdaImag: Vec[Double])

case class EigenDecompositionSymmetric(q: Mat[Double], lambdaReal: Vec[Double])

trait OpImpl {

  lazy val BLAS = com.github.fommil.netlib.BLAS.getInstance

  lazy val LAPACK = com.github.fommil.netlib.LAPACK.getInstance

  implicit def symmEigenValueTrunc =
    new MatUnaryOp1Scalar[EigValSymTrunc, Int, Vec[Double]] {
      def apply(m: Mat[Double], k: Int): Vec[Double] = {
        assert(m.numRows == m.numCols)
        val K = math.min(m.numRows, k)
        val a = m.contents.clone

        val vl = Array.ofDim[Double](1)
        val wr = Array.ofDim[Double](K)

        val workQuery = Array.ofDim[Double](1)
        val info = new org.netlib.util.intW(0)
        val outw = new org.netlib.util.intW(0)
        val ifail = Array.ofDim[Int](K)
        val iwork = Array.ofDim[Int](5 * m.numRows)

        LAPACK.dsyevx("N",
                      "I",
                      "U",
                      m.numRows,
                      a,
                      m.numRows,
                      0d,
                      0d,
                      1,
                      K,
                      0d,
                      outw,
                      wr,
                      vl,
                      m.numRows,
                      workQuery,
                      -1,
                      iwork,
                      ifail,
                      info)

        val work = Array.ofDim[Double](workQuery(0).toInt)

        LAPACK.dsyevx("N",
                      "I",
                      "U",
                      m.numRows,
                      a,
                      m.numRows,
                      0d,
                      0d,
                      m.numRows - K + 1,
                      m.numRows,
                      0d,
                      outw,
                      wr,
                      vl,
                      m.numRows,
                      work,
                      work.size,
                      iwork,
                      ifail,
                      info)

        val success = info.`val` == 0

        if (!success) throw new RuntimeException("Eigen decomposition failed")

        val wr2: Vec[Double] = (wr: Vec[Double]).reversed
        wr2
      }
    }

  implicit def symmEigenTrunc =
    new MatUnaryOp1Scalar[EigSTrunc, Int, EigenDecompositionSymmetric] {
      def apply(m: Mat[Double], k: Int): EigenDecompositionSymmetric = {
        assert(m.numRows == m.numCols)
        val K = math.min(m.numRows, k)
        val a = m.contents.clone

        val vl = Array.ofDim[Double](m.numRows * K)
        val wr = Array.ofDim[Double](K)

        val workQuery = Array.ofDim[Double](1)
        val info = new org.netlib.util.intW(0)
        val outw = new org.netlib.util.intW(0)
        val ifail = Array.ofDim[Int](K)
        val iwork = Array.ofDim[Int](5 * m.numRows)

        LAPACK.dsyevx("V",
                      "I",
                      "U",
                      m.numRows,
                      a,
                      m.numRows,
                      0d,
                      0d,
                      1,
                      K,
                      0d,
                      outw,
                      wr,
                      vl,
                      m.numRows,
                      workQuery,
                      -1,
                      iwork,
                      ifail,
                      info)

        val work = Array.ofDim[Double](workQuery(0).toInt)

        LAPACK.dsyevx("V",
                      "I",
                      "U",
                      m.numRows,
                      a,
                      m.numRows,
                      0d,
                      0d,
                      m.numRows - K + 1,
                      m.numRows,
                      0d,
                      outw,
                      wr,
                      vl,
                      m.numRows,
                      work,
                      work.size,
                      iwork,
                      ifail,
                      info)

        val success = info.`val` == 0

        if (!success) throw new RuntimeException("Eigen decomposition failed")

        val wr2: Vec[Double] = (wr: Vec[Double]).reversed
        val vl2 = new mat.MatDouble(K, m.numRows, vl)
          .takeRows((0 until K).reverse.toArray)

        EigenDecompositionSymmetric(vl2.T, wr2)
      }
    }

  implicit def symmEigen =
    new MatUnaryOp[EigS, EigenDecompositionSymmetric] {
      def apply(m: Mat[Double]): EigenDecompositionSymmetric = {
        assert(m.numRows == m.numCols)
        val a = m.contents.clone

        val wr = Array.ofDim[Double](m.numRows)

        val workQuery = Array.ofDim[Double](1)
        val info = new org.netlib.util.intW(0)

        LAPACK
          .dsyev("V", "U", m.numRows, a, m.numRows, wr, workQuery, -1, info)

        val work = Array.ofDim[Double](workQuery(0).toInt)

        LAPACK
          .dsyev("V", "U", m.numRows, a, m.numRows, wr, work, work.size, info)

        val success = info.`val` == 0

        if (!success) throw new RuntimeException("Eigen decomposition failed")

        val wr2: Vec[Double] = (wr: Vec[Double]).reversed
        val vl2 = new mat.MatDouble(m.numRows, m.numRows, a)
          .takeRows((0 until m.numRows).reverse.toArray)

        EigenDecompositionSymmetric(vl2.T, wr2)
      }
    }

  implicit def nonsymmEigen =
    new MatUnaryOp[EigNS, EigenDecompositionNonSymmetric] {
      def apply(m: Mat[Double]): EigenDecompositionNonSymmetric = {
        assert(m.numRows == m.numCols)
        val a = m.contents.clone

        val vl = Array.ofDim[Double](m.numRows * m.numRows)
        val wr = Array.ofDim[Double](m.numRows)
        val wi = Array.ofDim[Double](m.numRows)

        val workQuery = Array.ofDim[Double](1)
        val info = new org.netlib.util.intW(0)

        LAPACK.dgeev("V",
                     "N",
                     m.numRows,
                     a,
                     m.numRows,
                     wr,
                     wi,
                     vl,
                     m.numRows,
                     null,
                     1,
                     workQuery,
                     -1,
                     info)

        val work = Array.ofDim[Double](workQuery(0).toInt)

        LAPACK.dgeev("V",
                     "N",
                     m.numRows,
                     a,
                     m.numRows,
                     wr,
                     wi,
                     vl,
                     m.numRows,
                     null,
                     1,
                     work,
                     work.size,
                     info)

        val success = info.`val` == 0

        if (!success) throw new RuntimeException("Eigen decomposition failed")

        val reindex = array.argsort(wr).reverse
        val wr2: Vec[Double] = (wr: Vec[Double])(reindex)
        val wi2: Vec[Double] = (wi: Vec[Double])(reindex)
        val vl2 = new mat.MatDouble(m.numRows, m.numRows, vl).takeRows(reindex)

        EigenDecompositionNonSymmetric(vl2.T, wr2, wi2)
      }
    }

  // implicit def svdtrunc =
  //   new MatUnaryOp1Scalar[GeneralSVDTrunc, Int, SVDResult] {
  //     def apply(m: Mat[Double], k: Int): SVDResult = {
  //       val K = math.min(k, math.min(m.numRows, m.numCols))
  //
  //       /** Lapack gives us the SVD of the transpose
  //         * t(a) = v t(s) t(u)
  //         *   a  = u s t(v)
  //         */
  //       val cop = m.contents.clone
  //       val s = Array.ofDim[Double](K)
  //       val u = Array.ofDim[Double](m.numCols * K)
  //       val vt = Array.ofDim[Double](m.numRows * K)
  //       val lworkQuery = Array.ofDim[Double](1)
  //       val info = new org.netlib.util.intW(0)
  //       val ns = new org.netlib.util.intW(0)
  //
  //       // 1. Workspace query
  //       LAPACKE.dgesvdx(
  //         "A", // JOBU,
  //         "A", // JOBVT,
  //         "I", // RANGE
  //         m.numCols, // M,
  //         m.numRows, // N,
  //         cop, // A,
  //         m.numCols, // LDA,
  //         0d, //VL (not referenced)
  //         0d, //VU (not referenced)
  //         1, // IL,
  //         K, // IU
  //         ns, // NS
  //         s, // S,
  //         u, // U,
  //         m.numCols, // LDU,
  //         vt, // VT,
  //         m.numRows, // LDVT,
  //         lworkQuery, // WORK,
  //         -1, // LWORK,
  //         info // INFO
  //       )
  //
  //       val lwork = lworkQuery(0).toInt
  //       val work = Array.ofDim[Double](lwork)
  //
  //       LAPACKE.dgesvdx(
  //         "A", // JOBU,
  //         "A", // JOBVT,
  //         "I", // RANGE
  //         m.numCols, // M,
  //         m.numRows, // N,
  //         cop, // A,
  //         m.numCols, // LDA,
  //         0d, //VL (not referenced)
  //         0d, //VU (not referenced)
  //         1, // IL,
  //         K, // IU
  //         ns, // NS
  //         s, // S,
  //         u, // U,
  //         m.numCols, // LDU,
  //         vt, // VT,
  //         m.numRows, // LDVT,
  //         work, // WORK,
  //         lwork, // LWORK,
  //         info // INFO
  //       )
  //
  //       if (lapackInfoMethod.get(info) == 0) {
  //         val ut: Mat[Double] = new mat.MatDouble(m.numRows, K, vt)
  //         val v: Mat[Double] = new mat.MatDouble(K, m.numCols, u)
  //         val sigma: Vec[Double] = Vec(s: _*)
  //
  //         SVDResult(
  //           u = ut,
  //           sigma = sigma,
  //           vt = v
  //         )
  //       } else throw new RuntimeException("SVD Failed")
  //
  //     }
  //   }

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

      if (info.`val` == 0) {
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

  // private val lapackInfoMethod =
  //   java.lang.Class.forName("org.netlib.util.intW").getField("val")

  implicit def invertPD =
    new MatUnaryOp[InvertPDCholesky, Option[Mat[Double]]] {
      def apply(m: Mat[Double]): Option[Mat[Double]] = {

        val marray = m.contents
        val array = marray.clone
        val info = new org.netlib.util.intW(0)
        val info2 = new org.netlib.util.intW(0)

        LAPACK.dpotrf("L", m.numCols, array, m.numCols, info)
        LAPACK.dpotri("L", m.numCols, array, m.numCols, info2)

        if (info.`val` == 0 && info2.`val` == 0) {

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

          Some(new mat.MatDouble(m.numCols, m.numCols, array))

        } else if (info.`val` != 0) {
          if (info.`val` > 0) None
          else throw new DPotrfException(info.`val`)
        } else {
          throw new RuntimeException(
            "ERROR in dpotri info=" + info2.`val` + """
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

  implicit def trace =
    new MatUnaryOp[Trace, Double] {
      def apply(a: Mat[Double]): Double = {
        assert(a.numRows == a.numCols, "Trace of rectangular matrix")
        var s = 0.0
        var i = 0
        val d = a.contents
        while (i < a.numRows) {
          s += d(i * a.numRows + i)
          i += 1
        }
        s
      }
    }

  implicit def diag =
    new MatUnaryOp[Diag, Vec[Double]] {
      def apply(a: Mat[Double]): Vec[Double] = {
        val b = Array.ofDim[Double](math.min(a.numRows, a.numCols))
        var i = 0
        val d = a.contents
        while (i < b.size) {
          b(i) = d(i * a.numRows + i)
          i += 1
        }
        b
      }
    }

  implicit def ispd = new MatUnaryOp[TestPD, Boolean] {
    def apply(m: Mat[Double]): Boolean = {

      val marray = m.contents
      val array = marray.clone
      val info = new org.netlib.util.intW(0)
      LAPACK.dpotrf("L", m.numCols, array, m.numCols, info)
      info.`val` == 0
    }
  }

  implicit def svdtrunc =
    new MatUnaryOp1Scalar[GeneralSVDTrunc, Int, SVDResult] {
      def apply(m: Mat[Double], k: Int): SVDResult = {
        val K = math.min(k, math.min(m.numRows, m.numCols))

        if (m.numRows <= m.numCols) {
          val xxt = m.outerM
          val EigenDecompositionSymmetric(u, lambda) = xxt.eigSymm(K)
          val sigma = lambda.map(math.sqrt)
          val sigmainv = sigma.map(x => 1d / x)
          // inv(u) = t(u)
          val utm = u tmm m

          // inv(diag(sigma)) * utm
          val vt = Mat(utm.rows.zip(sigmainv.toSeq).map(x => x._1 * x._2): _*).T
          SVDResult(u, sigma, vt)
        } else {
          val xtx = m.innerM
          val EigenDecompositionSymmetric(v, lambda) = xtx.eigSymm(K)
          val sigma = lambda.map(math.sqrt)
          val sigmainv = sigma.map(x => 1d / x)

          val mv = m mm v
          // mv * inv(diag(sigma))
          val u = Mat(mv.cols.zip(sigmainv.toSeq).map(x => x._1 * x._2): _*)

          SVDResult(u, sigma, v.T)
        }

      }
    }

  implicit def singularValues =
    new MatUnaryOp1Scalar[SingularValues, Int, Vec[Double]] {
      def apply(m: Mat[Double], k: Int): Vec[Double] = {
        val K = math.min(k, math.min(m.numRows, m.numCols))

        if (m.numRows <= m.numCols) {
          val xxt = m.outerM
          xxt.eigenValuesSymm(K).map(math.sqrt)
        } else {
          val xtx = m.innerM
          xtx.eigenValuesSymm(K).map(math.sqrt)
        }

      }
    }

  /* D * M */
  implicit def multDiagFromLeft =
    new MatUnaryOp1Scalar[DiagxA, Vec[Double], Mat[Double]] {
      def apply(a: Mat[Double], b: Vec[Double]): Mat[Double] = {
        // assert(a.numRows == b.numCols,
        //        s"Incorrect dimensions ${a.numRows} ${b.numCols}")

        val result = Array.ofDim[Double](a.numCols * b.length)
        val ac = a.contents

        val I = b.length
        val J = a.numCols
        var i = 0
        var j = 0
        while (i < I) {
          while (j < J) {
            result(i * J + j) = a.contents(i * J + j) * b.raw(i)
            j += 1
          }
          j = 0
          i += 1
        }

        new mat.MatDouble(b.length, a.numCols, result)
      }
    }

  /* M * D */
  implicit def multDiagFromRight =
    new MatUnaryOp1Scalar[AxDiag, Vec[Double], Mat[Double]] {
      def apply(a: Mat[Double], b: Vec[Double]): Mat[Double] = {

        val result = Array.ofDim[Double](a.numRows * b.length)

        val ac = a.contents
        val I = a.numRows
        val J = b.length
        var i = 0
        var j = 0
        while (i < I) {
          while (j < J) {
            result(i * J + j) = a.contents(i * a.numCols + j) * b.raw(j)
            j += 1
          }
          j = 0
          i += 1
        }

        new mat.MatDouble(a.numRows, b.length, result)
      }
    }

}
