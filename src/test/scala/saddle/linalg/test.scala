package org.saddle.linalg

import org.saddle._
import org.scalatest.FunSuite

class VV1Suite extends FunSuite {
  test("1x3") {
    val m1 = Vec(1d, 2d, 3d)
    val m2 = Vec(4d, 5d, 6d)
    assert((m1 vv m2) == (4d + 10d + 18d))
  }

}

class MV1Suite extends FunSuite {
  test("1x3") {
    val m1 = Mat(Vec(1d, 2d, 3d), Vec(4d, 5d, 6d)).T
    val m2 = Vec(3d, 4d, 5d)
    assert(m1.mv(m2) == Vec(26d, 62d))
  }

}

class MVTSuite extends FunSuite {
  test("1x3") {
    val m1 = Mat(Vec(1d, 2d, 3d), Vec(4d, 5d, 6d))
    val m2 = Vec(3d, 4d, 5d)
    assert(m1.tmv(m2) == Vec(26d, 62d))
  }

}

class MMSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d))
    val m2 = Mat(Vec(3d, 4d))
    assert((m1 mm m2.T) == Mat(Vec(3d, 6d), Vec(4d, 8d)))
  }

  test("1x2 t") {
    val m1 = Mat(Vec(1d, 2d))
    val m2 = Mat(Vec(3d, 4d))
    assert((m1.T mm m2) == Mat(Vec(11d)))
  }

  test("2x3") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val r = m1 mm m1.T
    assert(r == Mat(Vec(35d, 44d), Vec(44d, 56d)))
  }

  test("2x3 t") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val r = m1.T mm m1
    assert(r == Mat(Vec(5d, 11d, 17d), Vec(11d, 25d, 39d), Vec(17d, 39d, 61d)))
  }
}

class MMCSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d))
    val m2 = Mat(Vec(3d, 4d))
    assert((m1 mmc (m2.T, mat.ones(2, 2))) == Mat(Vec(4d, 7d), Vec(5d, 9d)))
  }

}

class TMMCSuite extends FunSuite {
  test("1x2 ") {
    val m1 = Mat(Vec(1d, 2d)).T
    val m2 = Mat(Vec(3d, 4d))
    assert((m1.tmmc(m2.T, mat.ones(2, 2))) == Mat(Vec(4d, 7d), Vec(5d, 9d)))
  }
}

class TMMSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d)).T
    val m2 = Mat(Vec(3d, 4d))
    assert((m1 tmm m2.T) == Mat(Vec(3d, 6d), Vec(4d, 8d)))
  }

  test("1x2 t") {
    val m1 = Mat(Vec(1d, 2d)).T
    val m2 = Mat(Vec(3d, 4d))
    assert((m1.T tmm m2) == Mat(Vec(11d)))
  }

  test("2x3") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d)).T
    val r = m1 tmm m1
    assert(r == Mat(Vec(35d, 44d), Vec(44d, 56d)))
  }

  test("2x3 t") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val r = m1 tmm m1
    assert(r == Mat(Vec(5d, 11d, 17d), Vec(11d, 25d, 39d), Vec(17d, 39d, 61d)))
  }
}

class MMTCSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d))
    val m2 = Mat(Vec(3d, 4d))
    assert((m1.mmtc(m2, mat.ones(2, 2))) == Mat(Vec(4d, 7d), Vec(5d, 9d)))
  }
}

class InnerSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d))

    assert((m1.innerM) == Mat(Vec(5d)))
  }

  test("2x2") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d))

    assert((m1.innerM) == Mat(Vec(5d, 11d), Vec(11d, 25d)))
  }
}

class InnerSuitePlus extends FunSuite {

  test("2x2") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d)).innerMpC(1.0, 1.0, mat.ones(2, 2))

    assert((m1) == Mat(Vec(6d, 12d), Vec(12d, 26d)))
  }
}

class OuterSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d)).outerM

    assert(m1 == Mat(Vec(1d, 2d), Vec(2d, 4d)))
  }

}

class OuterSuitePlus extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d)).outerMpC(1.0, 1.0, Mat(Vec(1d, 1d), Vec(1d, 1d)))

    assert(m1 == Mat(Vec(2d, 3d), Vec(3d, 5d)))
  }

}

class MMTSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d))
    val m2 = Mat(Vec(3d, 4d))
    assert((m1 mmt m2) == Mat(Vec(3d, 6d), Vec(4d, 8d)))
  }

  test("1x2 t") {
    val m1 = Mat(Vec(1d, 2d))
    val m2 = Mat(Vec(3d, 4d)).T
    assert((m1.T mmt m2) == Mat(Vec(11d)))
  }

  test("2x3") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val r = m1 mmt m1
    assert(r == Mat(Vec(35d, 44d), Vec(44d, 56d)))
  }

  test("2x3 t") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val r = m1.T mmt m1.T
    assert(r == Mat(Vec(5d, 11d, 17d), Vec(11d, 25d, 39d), Vec(17d, 39d, 61d)))
  }
}

class TMMTCSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d)).T
    val m2 = Mat(Vec(3d, 4d))
    assert((m1.tmmtc(m2, mat.ones(2, 2))) == Mat(Vec(4d, 7d), Vec(5d, 9d)))
  }
}

class TMMTSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d)).T
    val m2 = Mat(Vec(3d, 4d))
    assert((m1 tmmt m2) == Mat(Vec(3d, 6d), Vec(4d, 8d)))
  }

  test("1x2 t") {
    val m1 = Mat(Vec(1d, 2d))
    val m2 = Mat(Vec(3d, 4d))
    assert((m1 tmmt m2.T) == Mat(Vec(11d)))
  }

  test("2x3") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val r = m1.T tmmt m1
    assert(r == Mat(Vec(35d, 44d), Vec(44d, 56d)))
  }

  test("2x3 t") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val r = m1 tmmt m1.T
    assert(r == Mat(Vec(5d, 11d, 17d), Vec(11d, 25d, 39d), Vec(17d, 39d, 61d)))
  }
}

class InvertGeneralSuite extends FunSuite {
  test("2x2") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d)).invert.roundTo(10)

    assert(m1 == Mat(Vec(-2d, 1d), Vec(1.5d, -0.5d)))
  }

}

class InvertPDSuite extends FunSuite {
  test("2x2") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d))
      .mmt(Mat(Vec(1d, 2d), Vec(3d, 4d)))
      .invertPD
      .get
      .roundTo(10)

    assert(Mat(Vec(1d, 2d), Vec(3d, 4d)).invertPD == None)

    assert(m1 == Mat(Vec(5d, -3.5), Vec(-3.5, 2.5d)))
  }

}

class SVD extends FunSuite {

  test("3x2") {
    val m = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val m1 = m.svd
    val sigma = Mat(mat.diag(m1.sigma).T.cols :+ Vec(0d, 0d): _*)
    val back = (m1.u mm sigma mm m1.vt)

    assert(back.roundTo(10) == m)
    assert(
      m1.u.roundTo(7) == Mat(Vec(-0.6196295, -0.7848945),
                             Vec(-0.7848945, 0.6196295)))

  }
}

class SVDTruncated extends FunSuite {

  test("3x2") {
    val m = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val mfull = m.svd(2)

    val sigma = Mat(mat.diag(mfull.sigma).T.cols: _*)
    val back = (mfull.u mm sigma mm mfull.vt)
    assert(back.roundTo(10) == m)
    val m1 = m.svd(1)

    assert(m1.u.roundTo(7) == Mat(Vec(0.6196295, 0.7848945)))

  }

  test("2x3") {
    val m = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d)).T
    val mfull = m.svd(2)

    val sigma = Mat(mat.diag(mfull.sigma).T.cols: _*)
    val back = (mfull.u mm sigma mm mfull.vt)
    assert(back.roundTo(10) == m)
    val m1 = m.svd(1)
    assert(m1.u.roundTo(4) == Mat(Vec(0.2298, 0.5247, 0.8196)))

  }
}

class SingulvarValues extends FunSuite {

  test("3x2") {
    val m = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val sv: Vec[Double] = m.singularValues(2)

    assert(sv == Vec(9.52551809156511, 0.5143005806586431))

  }

  test("2x3") {
    val m = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d)).T
    val sv = m.singularValues(2)

    assert(sv == Vec(9.52551809156511, 0.5143005806586431))

  }
}

class TraceSuite extends FunSuite {
  test("1x2") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d))

    assert((m1.trace) == 5d)
  }

}

class EigNSSuite extends FunSuite {
  test("2x2") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d))
    assert(
      m1.eigNonSymm.toString == EigenDecompositionNonSymmetric(
        Mat(Vec(-0.5658, -0.8246), Vec(-0.9094, 0.4160)),
        Vec(-0.3723, 5.3723).reversed,
        Vec(0d, 0d)).toString)
  }

}

class EigSSuite extends FunSuite {
  test("2x2") {
    val m1 = Mat(Vec(1d, 2d), Vec(2d, 1d))
    assert(
      m1.eigSymm.toString == EigenDecompositionSymmetric(
        Mat(Vec(0.7071, 0.7071), Vec(-0.7071, 0.7071)),
        Vec(3d, -1d)).toString)
  }

}

class EigSTruncSuite extends FunSuite {
  test("2x2") {
    val m1 = Mat(Vec(1d, 2d), Vec(2d, 1d))
    assert(
      m1.eigSymm(2).toString == EigenDecompositionSymmetric(
        Mat(Vec(0.7071, 0.7071), Vec(-0.7071, 0.7071)),
        Vec(3d, -1d)).toString)

    assert(
      m1.eigSymm(1).toString == EigenDecompositionSymmetric(
        Mat(Vec(0.7071, 0.7071)),
        Vec(3d)).toString)
  }

}

class DxMSuite extends FunSuite {
  test("2x2") {
    val m1 = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))

    assert(m1.mDiagFromLeft(Vec(10d)) == Mat(Vec(10d), Vec(30d), Vec(50d)))
    assert(
      m1.mDiagFromLeft(Vec(10d, 0.5)) == Mat(Vec(10d, 1d),
                                             Vec(30d, 2d),
                                             Vec(50d, 3d)))
  }

}

class MxDSuite extends FunSuite {
  test("2x2") {
    val m1 = Mat(Vec(1d, 2d, 3d), Vec(3d, 4d, 5d))

    assert(m1.mDiagFromRight(Vec(10d)) == Mat(Vec(10d, 20d, 30d)))
    assert(
      m1.mDiagFromRight(Vec(10d, 0.5)) == Mat(Vec(10d, 20d, 30d),
                                              Vec(1.5, 2d, 2.5)))
  }

}

class DiagAxAt extends FunSuite {
  test("2x3") {
    val x = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val diag = x.diagOuterM
    assert(diag.roundTo(3).col(0) == (x mm x.T).diag)
  }
  test("3x2") {
    val x = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d)).T
    val diag = x.diagOuterM
    assert(diag.roundTo(3).col(0) == (x mm x.T).diag)
  }

}

class DiagAtxA extends FunSuite {
  test("2x3") {
    val x = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val diag = x.diagInnerM
    assert(diag.roundTo(3).col(0) == (x.T mm x).diag)
  }
  test("3x2") {
    val x = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d)).T
    val diag = x.diagInnerM
    assert(diag.roundTo(3).col(0) == (x.T mm x).diag)
  }

}

class RowColSums extends FunSuite {
  test("2x3") {
    val x = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    assert(x.colSums.roundTo(3).col(0) == Vec(3d, 7d, 11d))
    assert(x.rowSums.roundTo(3).col(0) == Vec(9d, 12d))
  }
  test("3x2") {
    val x = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d)).T

    assert(x.rowSums.roundTo(3).col(0) == Vec(3d, 7d, 11d))
    assert(x.colSums.roundTo(3).col(0) == Vec(9d, 12d))
  }

}

class CholeskySuite extends FunSuite {
  test("2x2") {
    val a = Mat(Vec(1d, 2d), Vec(3d, 4d))
      .mmt(Mat(Vec(1d, 2d), Vec(3d, 4d)))
    assert(
      a.choleskyLower.get.roundTo(2) == Mat(Vec(3.16, 4.43), Vec(14d, 0.63)))

  }

}

class DeterminantPDSuite extends FunSuite {
  test("2x2") {
    val a = Mat(Vec(1d, 2d), Vec(3d, 4d))
      .mmt(Mat(Vec(1d, 2d), Vec(3d, 4d)))
    assert(math.abs(
      a.determinantPD.get - a.eigenValuesSymm(100).map(math.log10).sum) < 1E-10)

  }

}

class DiagXAInverseXtSuite extends FunSuite {
  test("2x3") {
    val x = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d)).T
    val a = Mat(Vec(1d, 2d), Vec(3d, 4d))
      .mmt(Mat(Vec(1d, 2d), Vec(3d, 4d)))

    val diag = a.diagInverseSandwich(x)
    assert(diag.get.roundTo(3) == Vec(1d, 1d, 5d))

    val cholesky = a.choleskyLower.get

    assert(
      cholesky.solveLowerTriangularForTransposed(x).get.T.roundTo(2) ==
        Mat(Vec(0.32, 0.95), Vec(0.95, -0.32), Vec(1.58, -1.58)))

    val z = cholesky.solveLowerTriangularForTransposed(x).get
    val diag2 = z.diagOuterM
    assert(diag.get.roundTo(3) == diag2.roundTo(3).col(0))

  }

}
class DiagXAInverseXtSuite2 extends FunSuite {
  test("2x3") {
    val x = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d)).T
    val a = Mat(Vec(1d, 2d), Vec(3d, 4d))
      .mmt(Mat(Vec(1d, 2d), Vec(3d, 4d)))

    val diag = a.diagInverseSandwich(x)
    assert(diag.get.roundTo(3) == Vec(1d, 1d, 5d))

    val cholesky = a.choleskyLower.get

    val choleskyT = {
      val a = cholesky.T.toArray
      a(2) = 0d
      Mat(2, 2, a)
    }

    assert(
      cholesky.T
        .solveUpperTriangularForTransposed((choleskyT mm a).T)
        .get
        .T
        .roundTo(2) == a)

  }

}
class SolveSuite extends FunSuite {
  test("2x2") {
    val b = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
    val a = Mat(Vec(1d, 2d), Vec(3d, 4d))
      .mmt(Mat(Vec(1d, 2d), Vec(3d, 4d)))
    val b1 = Mat(Vec(0.9999999999999964, 2d),
                 Vec(3.000000000000001, 4d),
                 Vec(5.000000000000007, 6d))

    val x = a.solve(b).get
    assert((a mm x) == b1)
    assert((a mm (a \ b).get) == b1)

  }

}
