package org.saddle

import org.saddle._

package object linalg extends OpImpl {
  implicit def pimp(m: Mat[Double]): Pimp = new Pimp(m)
}
