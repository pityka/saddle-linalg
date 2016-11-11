package org.saddle

import org.saddle._

package object linalg extends OpImpl {
  implicit def pimp(m: Mat[Double]): MatPimp = new MatPimp(m)
  implicit def pimp(m: Vec[Double]): VecPimp = new VecPimp(m)
}
