package org.saddle.linalg

import org.saddle._
import org.scalatest.FunSuite

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
      .roundTo(10)

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

// class SVDTruncated extends FunSuite {
//
//   test("3x2") {
//     val m = Mat(Vec(1d, 2d), Vec(3d, 4d), Vec(5d, 6d))
//     val m1 = m.svd(2)
//     // val sigma = Mat(mat.diag(m1.sigma).T.cols :+ Vec(0d, 0d): _*)
//     // val back = (m1.u mm sigma mm m1.vt)
//     println(m1)
//
//     println(m.svd)
//     // assert(back.roundTo(10) == m)
//     // assert(
//     //   m1.u.roundTo(7) == Mat(Vec(-0.6196295, -0.7848945),
//     //                          Vec(-0.7848945, 0.6196295)))
//
//   }
// }

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
